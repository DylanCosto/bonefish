//! RMQ-based anchor chaining
//!
//! An alternative to DP chaining that uses range-minimum queries on an
//! augmented treap to find optimal predecessors in O(log n) per anchor instead
//! of scanning backwards through a distance window. This is more efficient for
//! presets with large bandwidths such as `lr:hqae` and `asm`.
//!
//! The treap is keyed by query position and augmented with subtree-minimum
//! priority (negative DP score adjusted by position). An optional inner tree
//! with a tighter distance window (`rmq_inner_dist`) handles nearby anchors
//! that the RMQ might miss due to gap penalty interactions.
//!
//! Backtracking, chain extraction, and output format (score|count descriptors
//! plus reordered anchors) are shared with the DP chaining module via
//! [`chain_backtrack`].

use crate::align::sketch::Minimizer;
use crate::align::chain::{fast_log2, chain_backtrack};
use crate::align::sort::radix_sort_128x;
use crate::align::map::{ChainingParams, ChainingBuffers};

const NIL: u32 = u32::MAX;
/// Query-position bin width (log2) of the segment-tree rank lookup table.
const YBIN_SHIFT: usize = 6;

/// Arena-based AVL tree with subtree-min augmentation for O(log n) RMQ.
/// Replaces the previous treap with iterative insert/erase and compact nodes.
struct RmqTree {
    nodes: Vec<AvlNode>,
    root: u32,
    size: usize,
    free_head: u32, // free list for node reuse
    // Scratch buffers reused across insert_elem calls to avoid per-call allocation.
    insert_path: Vec<u32>,
    insert_stack: Vec<u8>,
}

#[derive(Clone, Copy)]
struct AvlNode {
    pri: f64,           // negative DP score (RMQ value to minimize)
    key_y: i32,         // query position (sort key)
    key_i: u32,         // anchor index (tiebreak key)
    left: u32,
    right: u32,
    parent: u32,
    sub_min_idx: u32,   // arena index of node with min pri in subtree
    balance: i8,        // AVL balance factor: -1, 0, +1
}

impl RmqTree {
    fn new() -> Self {
        RmqTree {
            nodes: Vec::with_capacity(256),
            root: NIL,
            size: 0,
            free_head: NIL,
            insert_path: Vec::with_capacity(64),
            insert_stack: Vec::with_capacity(64),
        }
    }

    #[inline]
    fn len(&self) -> usize { self.size }

    #[inline(always)]
    fn key_lt(ay: i32, ai: u32, by: i32, bi: u32) -> bool {
        ay < by || (ay == by && ai < bi)
    }

    fn alloc_node(&mut self, key_y: i32, key_i: u32, pri: f64) -> u32 {
        let idx;
        if self.free_head != NIL {
            idx = self.free_head;
            self.free_head = self.nodes[idx as usize].right; // free list uses right pointer
            self.nodes[idx as usize] = AvlNode {
                pri, key_y, key_i, left: NIL, right: NIL, parent: NIL,
                sub_min_idx: idx, balance: 0,
            };
        } else {
            idx = self.nodes.len() as u32;
            self.nodes.push(AvlNode {
                pri, key_y, key_i, left: NIL, right: NIL, parent: NIL,
                sub_min_idx: idx, balance: 0,
            });
        }
        idx
    }

    fn free_node(&mut self, idx: u32) {
        self.nodes[idx as usize].right = self.free_head;
        self.free_head = idx;
    }

    #[inline]
    fn update_sub_min(&mut self, idx: u32) {
        let i = idx as usize;
        let l = self.nodes[i].left;
        let r = self.nodes[i].right;
        let cur_pri = self.nodes[i].pri;
        let mut best = if l == NIL {
            idx
        } else {
            let l_min = self.nodes[l as usize].sub_min_idx;
            if cur_pri < self.nodes[l_min as usize].pri { idx } else { l_min }
        };
        if r != NIL {
            let r_min = self.nodes[r as usize].sub_min_idx;
            let best_pri = self.nodes[best as usize].pri;
            if best_pri >= self.nodes[r_min as usize].pri {
                best = r_min;
            }
        }
        self.nodes[i].sub_min_idx = best;
    }

    /// Set child[dir] of p to c, updating parent pointer.
    #[inline]
    fn set_child(&mut self, p: u32, dir: usize, c: u32) {
        if dir == 0 { self.nodes[p as usize].left = c; }
        else { self.nodes[p as usize].right = c; }
        if c != NIL { self.nodes[c as usize].parent = p; }
    }

    #[inline]
    fn child(&self, idx: u32, dir: usize) -> u32 {
        if dir == 0 { self.nodes[idx as usize].left } else { self.nodes[idx as usize].right }
    }

    /// Single rotation: rotate x up over its parent p.
    /// dir=0: left rotation (x is right child of p), dir=1: right rotation.
    ///
    /// `x` inherits p's old sub_min; p's sub_min is recomputed from
    /// (p.p[dir], q.p[dir]) — for dir=1 this is reversed from natural left/right,
    /// which matters for tie-breaking among equal-pri descendants.
    fn rotate(&mut self, p: u32, dir: usize) -> u32 {
        let x = self.child(p, 1 - dir);
        let mid = self.child(x, dir);
        let gp = self.nodes[p as usize].parent;
        let saved_p_sub = self.nodes[p as usize].sub_min_idx;
        let p_dir_child = self.child(p, dir);
        self.set_child(x, dir, p);
        self.set_child(p, 1 - dir, mid);
        self.nodes[x as usize].parent = gp;
        if gp != NIL {
            if self.nodes[gp as usize].left == p { self.nodes[gp as usize].left = x; }
            else { self.nodes[gp as usize].right = x; }
        }
        self.update_sub_min_with_children(p, p_dir_child, mid);
        self.nodes[x as usize].sub_min_idx = saved_p_sub;
        x
    }

    /// Fused double rotation. Layout:
    ///   before: (a, ((b, c)r, d)q)p  with dir=0
    ///   after:  ((a, b)p, (c, d)q)r
    /// `r` inherits p's old sub_min; p and q are recomputed from explicit
    /// post-rotation children using update_sub_min_with_children, since the
    /// arg order to that helper carries tie-break preference.
    fn rotate2(&mut self, p: u32, dir: usize) -> u32 {
        let opp = 1 - dir;
        let q = self.child(p, opp);
        let r = self.child(q, dir);
        let saved_p_sub = self.nodes[p as usize].sub_min_idx;
        let gp = self.nodes[p as usize].parent;
        let p_dir = self.child(p, dir);
        let r_dir = self.child(r, dir);
        let q_opp = self.child(q, opp);
        let r_opp = self.child(r, opp);
        self.update_sub_min_with_children(p, p_dir, r_dir);
        self.update_sub_min_with_children(q, q_opp, r_opp);
        self.nodes[r as usize].sub_min_idx = saved_p_sub;
        self.set_child(p, opp, r_dir);
        self.set_child(r, dir, p);
        self.set_child(q, dir, r_opp);
        self.set_child(r, opp, q);
        self.nodes[r as usize].parent = gp;
        if gp != NIL {
            if self.nodes[gp as usize].left == p { self.nodes[gp as usize].left = r; }
            else { self.nodes[gp as usize].right = r; }
        }
        r
    }

    /// Compute sub_min for `idx` from explicit children rather than reading
    /// `idx`'s current left/right. `q` is preferred over `r` on equal-pri.
    #[inline]
    fn update_sub_min_with_children(&mut self, idx: u32, q: u32, r: u32) {
        let cur_pri = self.nodes[idx as usize].pri;
        let mut best = idx;
        if q != NIL {
            let q_min = self.nodes[q as usize].sub_min_idx;
            if cur_pri >= self.nodes[q_min as usize].pri {
                best = q_min;
            }
        }
        if r != NIL {
            let r_min = self.nodes[r as usize].sub_min_idx;
            if self.nodes[best as usize].pri >= self.nodes[r_min as usize].pri {
                best = r_min;
            }
        }
        self.nodes[idx as usize].sub_min_idx = best;
    }

    /// Insert and rebalance. BST descent records the ancestor path and the
    /// deepest non-zero-balance ancestor `bp`. Then two independent walks:
    /// sub_min leaf->root with early break when the new node is no longer the
    /// subtree min, and balance bp->leaf. One rotation at bp if `|balance|>=2`.
    ///
    /// Sub_min and balance walks must not be interleaved: an interleaved walk
    /// that exits at the rotation point and then re-walks ancestors with
    /// natural-args `update_sub_min` overrides rotation-specific tie-break
    /// results at `new_root.parent` on equal-pri ties.
    fn insert_elem(&mut self, y: i32, i: usize, pri: f64) {
        let key_i = i as u32;
        let new_idx = self.alloc_node(y, key_i, pri);
        self.size += 1;

        if self.root == NIL {
            self.root = new_idx;
            return;
        }

        // BST descent. Record full ancestor path; bp_idx tracks the deepest
        // non-zero-balance ancestor seen so far. Balance walk later uses only
        // stack[bp_idx..]. Buffers are reused across calls to avoid per-call
        // allocation.
        let mut path = std::mem::take(&mut self.insert_path);
        let mut stack = std::mem::take(&mut self.insert_stack);
        path.clear();
        stack.clear();
        let mut bp_idx: usize = 0;
        let mut cur = self.root;
        loop {
            if self.nodes[cur as usize].balance != 0 {
                bp_idx = path.len();
            }
            let n = &self.nodes[cur as usize];
            let go_left = Self::key_lt(y, key_i, n.key_y, n.key_i);
            let dir = if go_left { 0u8 } else { 1u8 };
            let next = if go_left { n.left } else { n.right };
            path.push(cur);
            stack.push(dir);
            if next == NIL {
                self.set_child(cur, dir as usize, new_idx);
                break;
            }
            cur = next;
        }

        // Step 1: sub_min walk from leaf's parent up to root with early break.
        // alloc_node already set new leaf's sub_min_idx to itself.
        for k in (0..path.len()).rev() {
            let anc = path[k];
            self.update_sub_min(anc);
            if self.nodes[anc as usize].sub_min_idx != new_idx {
                break;
            }
        }

        // Step 2: balance walk from bp down to leaf's parent.
        for k in bp_idx..path.len() {
            let anc = path[k];
            let dir = stack[k];
            if dir == 0 {
                self.nodes[anc as usize].balance -= 1;
            } else {
                self.nodes[anc as usize].balance += 1;
            }
        }

        // Step 3: rotate at bp if balance hits +-2.
        let bp = path[bp_idx];
        let bp_balance = self.nodes[bp as usize].balance;
        if bp_balance >= 2 || bp_balance <= -2 {
            let heavy_dir = if bp_balance > 0 { 1usize } else { 0 };
            let heavy_child = self.child(bp, heavy_dir);
            let hb = self.nodes[heavy_child as usize].balance;
            let opposite = 1 - heavy_dir;
            let new_root;
            if (heavy_dir == 1 && hb >= 0) || (heavy_dir == 0 && hb <= 0) {
                new_root = self.rotate(bp, opposite);
                self.nodes[bp as usize].balance = if hb == 0 { if heavy_dir == 1 { 1 } else { -1 } } else { 0 };
                self.nodes[new_root as usize].balance = if hb == 0 { if heavy_dir == 1 { -1 } else { 1 } } else { 0 };
            } else {
                let inner = self.child(heavy_child, opposite);
                let inner_bal = self.nodes[inner as usize].balance;
                new_root = self.rotate2(bp, opposite);
                self.nodes[bp as usize].balance = if (heavy_dir == 1 && inner_bal > 0) || (heavy_dir == 0 && inner_bal < 0) {
                    if heavy_dir == 1 { -1 } else { 1 }
                } else { 0 };
                self.nodes[heavy_child as usize].balance = if (heavy_dir == 1 && inner_bal < 0) || (heavy_dir == 0 && inner_bal > 0) {
                    if heavy_dir == 1 { 1 } else { -1 }
                } else { 0 };
                self.nodes[new_root as usize].balance = 0;
                let _ = inner;
            }
            if self.root == bp { self.root = new_root; }
        }

        // Return scratch buffers to the tree for the next call.
        self.insert_path = path;
        self.insert_stack = stack;
    }

    fn erase(&mut self, y: i32, i: usize) -> bool {
        let key_i = i as u32;
        // Find the node
        let mut cur = self.root;
        while cur != NIL {
            let n = &self.nodes[cur as usize];
            if n.key_y == y && n.key_i == key_i { break; }
            cur = if Self::key_lt(y, key_i, n.key_y, n.key_i) { n.left } else { n.right };
        }
        if cur == NIL { return false; }

        self.size -= 1;

        // If two children, copy in-order successor's data and delete the successor instead
        let victim = cur;
        if self.nodes[cur as usize].left != NIL && self.nodes[cur as usize].right != NIL {
            let mut succ = self.nodes[cur as usize].right;
            while self.nodes[succ as usize].left != NIL {
                succ = self.nodes[succ as usize].left;
            }
            self.nodes[cur as usize].key_y = self.nodes[succ as usize].key_y;
            self.nodes[cur as usize].key_i = self.nodes[succ as usize].key_i;
            self.nodes[cur as usize].pri = self.nodes[succ as usize].pri;
            cur = succ;
        }

        // cur has at most one child. Record which side of parent it was on BEFORE splicing.
        let par = self.nodes[cur as usize].parent;
        let del_was_left = par != NIL && self.nodes[par as usize].left == cur;

        let child = if self.nodes[cur as usize].left != NIL {
            self.nodes[cur as usize].left
        } else {
            self.nodes[cur as usize].right
        };

        // Splice out cur
        if child != NIL { self.nodes[child as usize].parent = par; }
        if par == NIL {
            self.root = child;
        } else if del_was_left {
            self.nodes[par as usize].left = child;
        } else {
            self.nodes[par as usize].right = child;
        }

        self.free_node(cur);

        // Walk up from par rebalancing. Track which side shrunk.
        let mut p = par;
        let mut shrunk_left = del_was_left;
        while p != NIL {
            self.update_sub_min(p);
            let old_balance = self.nodes[p as usize].balance;
            // Adjust balance: left shrunk -> balance increases; right shrunk -> decreases
            let new_balance = if shrunk_left { old_balance + 1 } else { old_balance - 1 };
            self.nodes[p as usize].balance = new_balance;

            if new_balance == 1 || new_balance == -1 {
                // Height didn't change (was 0, now +-1), stop propagating
                break;
            } else if new_balance == 0 {
                // Height decreased, continue propagating up
                let pp = self.nodes[p as usize].parent;
                if pp != NIL {
                    shrunk_left = self.nodes[pp as usize].left == p;
                }
                p = pp;
            } else {
                // |new_balance| == 2, need rotation
                let heavy_dir = if new_balance > 0 { 1usize } else { 0 };
                let heavy_child = self.child(p, heavy_dir);
                let hb = self.nodes[heavy_child as usize].balance;
                let opposite = 1 - heavy_dir;

                if (heavy_dir == 1 && hb >= 0) || (heavy_dir == 0 && hb <= 0) {
                    // Single rotation
                    let new_root = self.rotate(p, opposite);
                    if hb == 0 {
                        // Height didn't change after rotation
                        self.nodes[p as usize].balance = if heavy_dir == 1 { 1 } else { -1 };
                        self.nodes[new_root as usize].balance = if heavy_dir == 1 { -1 } else { 1 };
                        if self.root == p { self.root = new_root; }
                        break; // height unchanged, stop
                    } else {
                        self.nodes[p as usize].balance = 0;
                        self.nodes[new_root as usize].balance = 0;
                        if self.root == p { self.root = new_root; }
                        // Height decreased, continue propagating
                        let pp = self.nodes[new_root as usize].parent;
                        if pp != NIL {
                            shrunk_left = self.nodes[pp as usize].left == new_root;
                        }
                        p = pp;
                    }
                } else {
                    // Double rotation (fused).
                    let inner = self.child(heavy_child, opposite);
                    let inner_bal = self.nodes[inner as usize].balance;
                    let new_root = self.rotate2(p, opposite);
                    self.nodes[p as usize].balance = if (heavy_dir == 1 && inner_bal > 0) || (heavy_dir == 0 && inner_bal < 0) {
                        if heavy_dir == 1 { -1 } else { 1 }
                    } else { 0 };
                    self.nodes[heavy_child as usize].balance = if (heavy_dir == 1 && inner_bal < 0) || (heavy_dir == 0 && inner_bal > 0) {
                        if heavy_dir == 1 { 1 } else { -1 }
                    } else { 0 };
                    self.nodes[new_root as usize].balance = 0;
                    if self.root == p { self.root = new_root; }
                    let _ = inner;
                    // Height decreased, continue propagating
                    let pp = self.nodes[new_root as usize].parent;
                    if pp != NIL {
                        shrunk_left = self.nodes[pp as usize].left == new_root;
                    }
                    p = pp;
                }
            }
        }
        // Continue updating sub_min for remaining ancestors
        while p != NIL {
            self.update_sub_min(p);
            p = self.nodes[p as usize].parent;
        }
        // Update sub_min from victim (whose key/pri was overwritten) up to root
        if victim != cur {
            let mut v = victim;
            while v != NIL {
                self.update_sub_min(v);
                v = self.nodes[v as usize].parent;
            }
        }
        true
    }

    /// Two-path LCA range minimum query: returns the element with minimum pri
    /// in the closed key interval [(lo_y, u32::MAX), (hi_y, 0)] — i.e. y > lo_y
    /// strictly, y < hi_y, or y == hi_y with i == 0.
    fn rmq(&self, lo_y: i32, hi_y: i32) -> Option<(i32, usize, f64)> {
        if self.root == NIL { return None; }

        // Key comparison: (y, i) ordering. lo = (lo_y, MAX), hi = (hi_y, 0).
        let lo_key = (lo_y, u32::MAX);
        let hi_key = (hi_y, 0u32);

        // Trace path from root to lo bound
        let mut path_lo: [(u32, i8); 64] = [(NIL, 0); 64]; // (node, cmp result)
        let mut len_lo = 0usize;
        let mut p = self.root;
        while p != NIL {
            let n = &self.nodes[p as usize];
            let nk = (n.key_y, n.key_i);
            let cmp = if lo_key < nk { -1i8 } else if lo_key > nk { 1 } else { 0 };
            path_lo[len_lo] = (p, cmp);
            len_lo += 1;
            if cmp < 0 { p = n.left; }
            else if cmp > 0 { p = n.right; }
            else { break; }
        }

        // Trace path from root to hi bound
        let mut path_hi: [(u32, i8); 64] = [(NIL, 0); 64];
        let mut len_hi = 0usize;
        p = self.root;
        while p != NIL {
            let n = &self.nodes[p as usize];
            let nk = (n.key_y, n.key_i);
            let cmp = if hi_key < nk { -1i8 } else if hi_key > nk { 1 } else { 0 };
            path_hi[len_hi] = (p, cmp);
            len_hi += 1;
            if cmp < 0 { p = n.left; }
            else if cmp > 0 { p = n.right; }
            else { break; }
        }

        // Find LCA: first shared node where lo goes left/equal and hi goes right/equal
        let mut lca = 0usize;
        let mut found_lca = false;
        for i in 0..len_lo.min(len_hi) {
            if path_lo[i].0 == path_hi[i].0 && path_lo[i].1 <= 0 && path_hi[i].1 >= 0 {
                lca = i;
                found_lca = true;
                break;
            }
            if path_lo[i].0 != path_hi[i].0 { break; }
        }
        if !found_lca { return None; }

        let mut min_idx = path_lo[lca].0; // start with LCA node

        // Scan lo-path below LCA: nodes where we went left/equal are in range.
        // Their right subtrees are fully in range.
        for (node, cmp) in path_lo.iter().take(len_lo).skip(lca + 1) {
            if *cmp <= 0 {
                // This node is in range (>= lo)
                if self.nodes[*node as usize].pri < self.nodes[min_idx as usize].pri {
                    min_idx = *node;
                }
                // Its right subtree is fully in range
                let rc = self.nodes[*node as usize].right;
                if rc != NIL {
                    let rc_min = self.nodes[rc as usize].sub_min_idx;
                    if self.nodes[rc_min as usize].pri < self.nodes[min_idx as usize].pri {
                        min_idx = rc_min;
                    }
                }
            }
        }

        // Scan hi-path below LCA: nodes where we went right/equal are in range.
        // Their left subtrees are fully in range.
        for (node, cmp) in path_hi.iter().take(len_hi).skip(lca + 1) {
            if *cmp >= 0 {
                if self.nodes[*node as usize].pri < self.nodes[min_idx as usize].pri {
                    min_idx = *node;
                }
                let lc = self.nodes[*node as usize].left;
                if lc != NIL {
                    let lc_min = self.nodes[lc as usize].sub_min_idx;
                    if self.nodes[lc_min as usize].pri < self.nodes[min_idx as usize].pri {
                        min_idx = lc_min;
                    }
                }
            }
        }

        let n = &self.nodes[min_idx as usize];
        Some((n.key_y, n.key_i as usize, n.pri))
    }

}

/// Segment-tree RMQ over the query-position ranks of all anchors of one chaining call. Holds the
/// same (key, priority) set as the AVL tree and answers the same range-minimum queries; only the
/// choice among exactly equal priorities can differ.
struct SegRmq {
    ys: Vec<i32>,       // query position by rank
    ybin: Vec<u32>,     // ybin[k] = first rank with query position >= k << YBIN_SHIFT
    is: Vec<u32>,       // anchor index by rank
    rank: Vec<u32>,     // rank by anchor index
    pri: Vec<f64>,      // 2 * cap, 1-based heap layout
    idx: Vec<u32>,
    cap: usize,
    size: usize,
}

impl SegRmq {
    fn new(a: &[Minimizer]) -> Self {
        let n = a.len();
        let mut keys: Vec<u64> = (0..n).map(|j| ((a[j].query_pos() as u32 as u64) << 32) | j as u64).collect();
        crate::align::sort::radix_sort_u64(&mut keys); // distinct keys (anchor index in the low bits)
        let mut ys = Vec::with_capacity(n);
        let mut is = Vec::with_capacity(n);
        let mut rank = vec![0u32; n];
        for (r, &k) in keys.iter().enumerate() {
            ys.push((k >> 32) as u32 as i32);
            is.push(k as u32);
            rank[k as u32 as usize] = r as u32;
        }
        let cap = n.next_power_of_two().max(1);
        let ymax = ys.last().copied().unwrap_or(0).max(0) as usize;
        let mut ybin = Vec::with_capacity((ymax >> YBIN_SHIFT) + 2);
        let mut r = 0usize;
        // (No table if a query position were negative: `bound` then searches all ranks.)
        let nbins = if ys.first().is_some_and(|&y| y < 0) { 0 } else { (ymax >> YBIN_SHIFT) + 2 };
        for k in 0..nbins {
            let t = (k << YBIN_SHIFT) as i64;
            while r < ys.len() && (ys[r] as i64) < t { r += 1; }
            ybin.push(r as u32);
        }
        SegRmq { ys, ybin, is, rank, pri: vec![f64::INFINITY; 2 * cap], idx: vec![NIL; 2 * cap], cap, size: 0 }
    }

    #[inline]
    fn len(&self) -> usize { self.size }

    #[inline(always)]
    fn pull(&mut self, mut p: usize) {
        while p > 1 {
            p >>= 1;
            // Child index by arithmetic, not a (mispredicting) branch.
            let c = 2 * p + (self.pri[2 * p + 1] < self.pri[2 * p]) as usize;
            let (np, ni) = (self.pri[c], self.idx[c]);
            // Unchanged node: its ancestors (functions of their children only) are unchanged too.
            if np.to_bits() == self.pri[p].to_bits() && ni == self.idx[p] { break; }
            self.pri[p] = np;
            self.idx[p] = ni;
        }
    }

    /// First rank whose query position is > y (`le`) or >= y (`!le`): the coarse bin table
    /// narrows it to the ranks of one 64-bp bin, then a short binary search finishes.
    #[inline(always)]
    fn bound(&self, y: i32, le: bool) -> usize {
        let ys = &self.ys;
        let t = y as i64 + le as i64; // first rank with query position >= t
        let (mut base, end) = if self.ybin.is_empty() {
            (0, ys.len())
        } else {
            if t <= 0 { return 0; }
            let k = (t >> YBIN_SHIFT) as usize;
            if k + 1 >= self.ybin.len() { return ys.len(); }
            (self.ybin[k] as usize, self.ybin[k + 1] as usize)
        };
        if t > i32::MAX as i64 { return ys.len(); }
        let t = t as i32;
        let mut n = end - base;
        while n > 0 {
            let half = n / 2;
            if ys[base + half] < t { base += half + 1; n -= half + 1; } else { n = half; }
        }
        base
    }

    #[inline]
    fn insert_elem(&mut self, j: usize, pri: f64) {
        let p = self.cap + self.rank[j] as usize;
        self.pri[p] = pri;
        self.idx[p] = j as u32;
        self.size += 1;
        self.pull(p);
    }

    #[inline]
    fn erase(&mut self, j: usize) {
        let p = self.cap + self.rank[j] as usize;
        if self.idx[p] == NIL { return; }
        self.pri[p] = f64::INFINITY;
        self.idx[p] = NIL;
        self.size -= 1;
        self.pull(p);
    }

    /// Minimum over keys (y, i) with lo_y < y < hi_y, or y == hi_y and i == 0 (the AVL query
    /// bounds (lo_y, MAX) .. (hi_y, 0)).
    #[inline]
    fn rmq(&self, lo_y: i32, hi_y: i32) -> Option<(i32, usize, f64)> {
        let mut l = self.bound(lo_y, true);
        let mut r = self.bound(hi_y, false);
        if r < self.ys.len() && self.ys[r] == hi_y && self.is[r] == 0 { r += 1; }
        if l >= r { return None; }
        let (mut best_p, mut best_i) = (f64::INFINITY, NIL);
        l += self.cap;
        r += self.cap;
        // Branch-free over the parity bits (which side takes a node is data-dependent);
        // a node is read even when it is not taken, but only counts when taken. The left
        // node is tested before the right one, as in the branchy form.
        while l < r {
            let lt = l & 1 == 1;
            let better = lt & (self.pri[l] < best_p);
            best_i = if better { self.idx[l] } else { best_i };
            best_p = f64::from_bits(if better { self.pri[l].to_bits() } else { best_p.to_bits() });
            l += lt as usize;
            let rt = r & 1 == 1;
            r -= rt as usize;
            let better = rt & (self.pri[r] < best_p);
            best_i = if better { self.idx[r] } else { best_i };
            best_p = f64::from_bits(if better { self.pri[r].to_bits() } else { best_p.to_bits() });
            l >>= 1;
            r >>= 1;
        }
        if best_i == NIL { return None; }
        let j = best_i as usize;
        Some((self.ys[self.rank[j] as usize], j, best_p))
    }
}

pub static RMQ_CHECK_DIFF: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static RMQ_CHECK_DIFF_PRI: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
pub static RMQ_CHECK_N: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// The inner (short-range) RMQ "tree" is only ever iterated in reverse key order, never range-
/// minimum queried, so a sorted vector of (y, i) keys gives the identical iteration order with
/// cheaper insert/erase (it holds the anchors of a short reference window).
struct SortedKeys {
    keys: Vec<u64>,
}

impl SortedKeys {
    fn new() -> Self { SortedKeys { keys: Vec::with_capacity(256) } }

    #[inline(always)]
    fn key(y: i32, i: usize) -> u64 { ((y as u32 as u64) << 32) | i as u32 as u64 }

    #[inline]
    fn len(&self) -> usize { self.keys.len() }

    #[inline]
    fn insert_elem(&mut self, y: i32, i: usize) {
        let k = Self::key(y, i);
        let pos = self.keys.partition_point(|&x| x < k);
        self.keys.insert(pos, k);
    }

    #[inline]
    fn erase(&mut self, y: i32, i: usize) {
        let k = Self::key(y, i);
        if let Ok(pos) = self.keys.binary_search(&k) {
            self.keys.remove(pos);
        }
    }

    /// Elements with key_y <= start_y in decreasing (y, i) order, as (y, i).
    #[inline]
    fn iter_rev_le(&self, start_y: i32) -> impl Iterator<Item = (i32, usize)> + '_ {
        let end = self.keys.partition_point(|&x| (x >> 32) as u32 as i32 <= start_y);
        self.keys[..end].iter().rev().map(|&k| ((k >> 32) as u32 as i32, (k & 0xffff_ffff) as usize))
    }
}

/// Simplified score computation for RMQ chaining.
/// Returns (score, is_exact, width). `width` is the diagonal deviation
/// `|ref_diff - query_diff|`, used to gate by bandwidth in the chaining loop.
#[inline(always)]
fn compute_chain_score_simple(
    ai: &Minimizer,
    aj: &Minimizer,
    chn_pen_gap: f32,
    chn_pen_skip: f32,
) -> (i32, bool, i32) {
    let query_diff = ai.query_pos().wrapping_sub(aj.query_pos());
    let ref_diff = ai.ref_pos().wrapping_sub(aj.ref_pos());
    let gap_width = if ref_diff > query_diff { ref_diff - query_diff } else { query_diff - ref_diff };
    let min_diff = if ref_diff < query_diff { ref_diff } else { query_diff };
    let q_span = aj.query_span();

    let mut sc = if q_span < min_diff { q_span } else { min_diff };
    let exact = gap_width == 0 && min_diff <= q_span;

    if gap_width > 0 || query_diff > q_span {
        let lin_pen = chn_pen_gap * (gap_width as f32) + chn_pen_skip * (min_diff as f32);
        let log_pen = if gap_width >= 1 { fast_log2((gap_width + 1) as f32) } else { 0.0f32 };
        sc -= (lin_pen + 0.5f32 * log_pen) as i32;
    }

    (sc, exact, gap_width)
}

/// RMQ-based chaining: uses a range-max query tree over score-by-position to
/// pick the best predecessor in O(log n) per anchor instead of the DP path's
/// O(n) inner scan. Wins on large bandwidths / long-read chains.
///
/// Arguments:
/// - opt: chaining parameters (max_gap, rmq_inner_dist, bandwidth, max_chain_skip, etc.)
/// - a: input anchors (modified in place)
/// - ctx: reusable map context (provides DP buffers)
///
/// Returns: (u, a_new) where u is chain scores/counts and a_new is the reordered anchors
pub fn chain_anchors_rmq(
    opt: &ChainingParams,
    a: &mut [Minimizer],
    ctx: &mut ChainingBuffers,
) -> (Vec<u64>, Vec<Minimizer>) {
    let n = a.len();
    if n == 0 {
        return (Vec::new(), Vec::new());
    }

    let bw = opt.bandwidth;
    let max_dist = if opt.max_gap < bw { bw } else { opt.max_gap };
    let max_dist_inner = if opt.rmq_inner_dist < 0 { 0 }
        else if opt.rmq_inner_dist > max_dist { max_dist }
        else { opt.rmq_inner_dist };
    let max_drop = bw;

    let mut predecessors = std::mem::take(&mut ctx.predecessors);
    let mut bt_candidates = std::mem::take(&mut ctx.bt_candidates);
    let mut scores = std::mem::take(&mut ctx.scores);
    let mut peak_scores = std::mem::take(&mut ctx.peak_scores);
    let mut visited = std::mem::take(&mut ctx.visited);
    predecessors.resize(n, 0i32);
    scores.resize(n, 0i32);
    peak_scores.resize(n, 0i32);
    // visited uses sentinel comparison (visited[j] == i) so must be zeroed
    visited.clear(); visited.resize(n, 0i32);

    let seg = opt.seg_rmq || *crate::align::env_flags::EXP_SEG_RMQ;
    let check = *crate::align::env_flags::EXP_RMQ_CHECK;
    let use_avl = !seg || check;
    let mut root = RmqTree::new();
    let mut seg_root = if seg || check { Some(SegRmq::new(a)) } else { None };
    let mut root_inner = if max_dist_inner > 0 { Some(SortedKeys::new()) } else { None };

    let mut window_start: usize = 0;
    let mut inner_window_start: usize = 0;
    let mut insert_from: usize = 0;
    let mut _mmax_f = 0;

    for i in 0..n {
        let mut best_predecessor: i64 = -1;
        let q_span = a[i].query_span();
        let mut best_score = q_span;

        // Add in-range anchors (when position changes)
    
        if insert_from < i && a[insert_from].x != a[i].x {
            for j in insert_from..i {
                let pri = -(scores[j] as f64 + 0.5 * opt.chn_pen_gap as f64 * ((a[j].ref_pos() + a[j].query_pos()) as f64));
                if use_avl { root.insert_elem(a[j].query_pos(), j, pri); }
                if let Some(sr) = &mut seg_root { sr.insert_elem(j, pri); }
                if let Some(inner) = &mut root_inner {
                    inner.insert_elem(a[j].query_pos(), j);
                }
            }
            insert_from = i;
        }

        // Remove anchors out of range from root
        while window_start < i && (
            a[i].ref_id_strand() != a[window_start].ref_id_strand() ||
            a[i].x > a[window_start].x + max_dist as u64 ||
            (if use_avl { root.len() } else { seg_root.as_ref().map_or(0, |sr| sr.len()) }) > opt.rmq_size_cap as usize
        ) {
            if use_avl { root.erase(a[window_start].query_pos(), window_start); }
            if let Some(sr) = &mut seg_root { sr.erase(window_start); }
            window_start += 1;
        }

        // Remove anchors out of range from root_inner
        if let Some(inner) = &mut root_inner {
            while inner_window_start < i && (
                a[i].ref_id_strand() != a[inner_window_start].ref_id_strand() ||
                a[i].x > a[inner_window_start].x + max_dist_inner as u64 ||
                inner.len() > opt.rmq_size_cap as usize
            ) {
                inner.erase(a[inner_window_start].query_pos(), inner_window_start);
                inner_window_start += 1;
            }
        }

        // RMQ query on root tree
        let lo_y = a[i].query_pos() - max_dist;
        let hi_y = a[i].query_pos();

        let hit = if use_avl {
            let h = root.rmq(lo_y, hi_y);
            if check {
                let h2 = seg_root.as_ref().unwrap().rmq(lo_y, hi_y);
                RMQ_CHECK_N.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if h.map(|x| x.1) != h2.map(|x| x.1) {
                    RMQ_CHECK_DIFF.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if h.map(|x| x.2) != h2.map(|x| x.2) { RMQ_CHECK_DIFF_PRI.fetch_add(1, std::sync::atomic::Ordering::Relaxed); }
                }
            }
            h
        } else {
            seg_root.as_ref().unwrap().rmq(lo_y, hi_y)
        };
        if let Some((_q_y, q_i, _q_pri)) = hit {
            let j = q_i;
            let (sc, exact, width) = compute_chain_score_simple(&a[i], &a[j], opt.chn_pen_gap, opt.chn_pen_skip);
            let total_sc = scores[j] + sc;
            if width <= bw && total_sc > best_score {
                best_score = total_sc;
                best_predecessor = j as i64;
            }

            // If not exact match, also search inner tree for close matches
            if !exact && let Some(ref inner) = root_inner && a[i].query_pos() > 0 {
                let mut skip_count = 0;
                for (q_y, q_i) in inner.iter_rev_le(a[i].query_pos() - 1) {
                    if q_y < a[i].query_pos() - max_dist_inner {
                        break;
                    }
                    let j = q_i;
                    let (sc, _, width) = compute_chain_score_simple(&a[i], &a[j], opt.chn_pen_gap, opt.chn_pen_skip);
                    if width <= bw {
                        let total_sc = scores[j] + sc;
                        if total_sc > best_score {
                            best_score = total_sc;
                            best_predecessor = j as i64;
                            if skip_count > 0 { skip_count -= 1; }
                        } else if visited[j] == i as i32 {
                            skip_count += 1;
                            if skip_count > opt.max_chain_skip {
                                break;
                            }
                        }
                        if predecessors[j] >= 0 {
                            visited[predecessors[j] as usize] = i as i32;
                        }
                    }
                }
            }
        }

        // Set max
        scores[i] = best_score;
        predecessors[i] = best_predecessor as i32;
        peak_scores[i] = if best_predecessor >= 0 && peak_scores[best_predecessor as usize] > best_score { peak_scores[best_predecessor as usize] } else { best_score };
        if _mmax_f < best_score { _mmax_f = best_score; }
    }

    // Backtrack to extract chains
    let (u, n_u, n_v) = chain_backtrack(n, &scores, &predecessors, &mut peak_scores, &mut visited, &mut bt_candidates, opt.min_cnt, opt.min_chain_score, max_drop);

    if n_u == 0 {
        ctx.predecessors = predecessors; ctx.scores = scores; ctx.peak_scores = peak_scores; ctx.visited = visited; ctx.bt_candidates = bt_candidates;
        return (Vec::new(), Vec::new());
    }

    // compact_a logic: reorder anchors into per-chain contiguous runs.
    // Step 1: Write chain anchors to b[] in forward order
    let mut b: Vec<Minimizer> = Vec::with_capacity(n_v);
    let mut k = 0;
    for &u_val in &u[..n_u] {
        let ni = (u_val & 0xFFFFFFFF) as usize;
        let k0 = k;
        for j in 0..ni {
            let idx = peak_scores[k0 + (ni - j - 1)] as usize;
            b.push(a[idx]);
            k += 1;
        }
    }

    // Step 2: Sort chains by target position of their first anchor.
    let mut w: Vec<Minimizer> = Vec::with_capacity(n_u);
    let mut k_pos = 0usize;
    for (i, &u_val) in u[..n_u].iter().enumerate() {
        let ni = (u_val & 0xFFFFFFFF) as usize;
        w.push(Minimizer {
            x: b[k_pos].x,
            y: ((k_pos as u64) << 32) | (i as u64),
        });
        k_pos += ni;
    }
    radix_sort_128x(&mut w);

    // Step 3: Reorder u[] and anchors according to sorted order
    let mut u2: Vec<u64> = Vec::with_capacity(n_u);
    let mut b2: Vec<Minimizer> = Vec::with_capacity(n_v);
    for &w_val in &w[..n_u] {
        let j = (w_val.y & 0xFFFFFFFF) as usize; // original chain index
        let offset = (w_val.y >> 32) as usize;    // offset in b[]
        let ni = (u[j] & 0xFFFFFFFF) as usize;
        u2.push(u[j]);
        for idx in 0..ni {
            b2.push(b[offset + idx]);
        }
    }

    ctx.predecessors = predecessors; ctx.scores = scores; ctx.peak_scores = peak_scores; ctx.visited = visited; ctx.bt_candidates = bt_candidates;
    (u2, b2)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rmq_tree_basic() {
        let mut tree = RmqTree::new();

        // Insert some elements
        tree.insert_elem(100, 0, -10.0);
        tree.insert_elem(200, 1, -20.0);
        tree.insert_elem(150, 2, -15.0);

        assert_eq!(tree.len(), 3);

        // RMQ closed interval [(lo_y, INT32_MAX), (hi_y, 0)]:
        //   y > lo_y (exclusive lower), y < hi_y (all), y == hi_y only if i == 0
        // Range (99, 201) includes y=100,150,200 → should return i=1 (pri=-20.0)
        let result = tree.rmq(99, 201);
        assert!(result.is_some());
        assert_eq!(result.unwrap().1, 1); // pri=-20.0 is lowest

        // Range (99, 150) includes only y=100 → should return i=0
        // Element at (150, i=2) excluded because i=2 > 0
        let result = tree.rmq(99, 150);
        assert!(result.is_some());
        assert_eq!(result.unwrap().1, 0);

        // Range (100, 200) excludes y=100 (lower bound exclusive), includes y=150 → i=2
        // Element at (200, i=1) excluded because i=1 > 0
        let result = tree.rmq(100, 200);
        assert!(result.is_some());
        assert_eq!(result.unwrap().1, 2);

        // Erase and check
        tree.erase(150, 2);
        assert_eq!(tree.len(), 2);
    }

    #[test]
    fn test_chain_rmq_simple() {
        let k = 15;
        let span_mask = (k as u64) << 32;
        let mut anchors: Vec<Minimizer> = Vec::new();

        // Linear chain of 3 anchors
        anchors.push(Minimizer { x: 100, y: span_mask | 100 });
        anchors.push(Minimizer { x: 120, y: span_mask | 120 });
        anchors.push(Minimizer { x: 150, y: span_mask | 150 });

        let opt = ChainingParams {
            min_cnt: 1,
            min_chain_score: 10,
            max_gap: 500,
            max_gap_ref: -1,
            max_dist_x: 500,
            max_dist_y: 500,
            bandwidth: 500,
            bandwidth_long: 500,
            max_chain_skip: 25,
            max_chain_iter: 5000,
            chn_pen_gap: 0.5,
            chn_pen_skip: 0.5,
            chain_gap_scale: 0.8,
            rmq_rescue_size: 1000,
            rmq_rescue_ratio: 0.1,
            rmq_inner_dist: 100,
            rmq_size_cap: 500,
            seg_rmq: false, rescue_join: false, rescue_max_density: 0,
        };
        let mut bufs = ChainingBuffers::new();
        let (u, _chains) = chain_anchors_rmq(
            &opt, &mut anchors, &mut bufs,
        );

        // Should produce at least one chain
        assert!(!u.is_empty());
    }
}
