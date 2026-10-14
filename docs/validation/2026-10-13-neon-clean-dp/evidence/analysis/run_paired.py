import datetime
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path('${SESSION}')
OUT = ROOT / 'outputs/bonefish-neon-clean-dp'
CONFIG = json.loads((OUT / 'EXPERIMENT_METADATA.json').read_text())
SCRATCH = Path(CONFIG['scratch'])
SCRATCH.mkdir(parents=True, exist_ok=True)
RUNS = []
HASHES = {}

def status(stage, **fields):
    value = dict(stage=stage, updated_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                 completed_runs=len(RUNS), total_runs=len(CONFIG['order']), **fields)
    (OUT / 'STATUS.json').write_text(json.dumps(value, indent=2) + '\n')
    print(json.dumps(value), flush=True)

try:
    for variant, trial in CONFIG['order']:
        name = f'{variant}-{trial}'
        if shutil.disk_usage(SCRATCH).free < 12 * (1 << 30):
            raise RuntimeError('Less than 12 GiB scratch space free before trial')
        status('running', run=name)
        command = [sys.executable, '-B', CONFIG['runner'], CONFIG[variant], '--reference', CONFIG['reference'],
                   '--reads', CONFIG['reads'], '--threads', str(CONFIG['threads']), '--variant', variant,
                   '--trial', str(trial), '--output-dir', str(SCRATCH)]
        if variant == 'control':
            command.extend(['--expected-sam-sha256', CONFIG['original_control_sam_sha256']])
        with (OUT / f'{name}.driver.log').open('wb') as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        for extension in ('.json', '.log'):
            path = SCRATCH / (name + extension)
            if path.exists():
                shutil.copy2(path, OUT / path.name)
        if result.returncode:
            raise RuntimeError(f'{name} failed with exit {result.returncode}; SAM retained')
        run = json.loads((SCRATCH / f'{name}.json').read_text())
        assert run['all_checks_passed'] and run['binary_sha256'] == CONFIG[variant + '_sha256']
        assert not run['warnings']
        digest = run['sam_sha256_without_pg']
        if variant in HASHES and digest != HASHES[variant]:
            raise RuntimeError(f'{variant} SAM changed between trials; SAM retained')
        HASHES[variant] = digest
        RUNS.append(run)
        # These are owned scratch outputs. Keep the first complete pair for detailed
        # comparison; later SAMs can be discarded only after full validation and receipt.
        if trial != 1:
            (SCRATCH / f'{name}.sam').unlink()
        status('trial_complete', run=name, seconds=run['elapsed_seconds_including_index_and_sam'],
               sam_sha256=digest, free_scratch_bytes=shutil.disk_usage(SCRATCH).free)
    (OUT / 'RUNS.json').write_text(json.dumps(RUNS, indent=2) + '\n')
    status('complete', variant_sam_sha256=HASHES)
except BaseException as error:
    status('failed', error=str(error))
    raise
