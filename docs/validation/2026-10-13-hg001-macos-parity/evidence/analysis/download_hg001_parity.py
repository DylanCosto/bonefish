from pathlib import Path
import datetime
import hashlib
import json
import subprocess
import time

SESSION = Path('${SESSION}')
OUT = SESSION / 'outputs/bonefish-hg001-parity'
WORK = SESSION / 'work/bonefish-hg001-parity'
DEST = WORK / 'data/hg001_wg_1x.fq.gz'
REMOTE = 'gdrive:hg001_wg_1x.fq.gz'
FILE_ID = '1L2O5CAHXbbypQxm52ZLbp6xAqCYYS7XH'
OUT.mkdir(parents=True, exist_ok=True)
DEST.parent.mkdir(parents=True, exist_ok=True)

def metadata():
    result = subprocess.run(['rclone', 'lsjson', REMOTE, '--stat', '--hash'],
                            text=True, capture_output=True, check=True)
    row = json.loads(result.stdout)
    assert row['ID'] == FILE_ID and row['Size'] == 1945506828
    assert row['Hashes']['sha256'] == '2edff9a3834d5d55963da198bad3717064e0ed5232a037d6783c84fbc8ffb288'
    return row

started = datetime.datetime.now(datetime.timezone.utc).isoformat()
before = metadata()
(OUT / 'DRIVE_SOURCE.json').write_text(json.dumps(before, indent=2) + '\n')
command = ['rclone', 'copyto', REMOTE, str(DEST), '--checksum', '--immutable',
           '--stats', '15s', '--stats-one-line', '--log-level', 'INFO',
           '--log-file', str(WORK / 'download.log')]
status_path = OUT / 'DOWNLOAD_STATUS.json'
status_path.write_text(json.dumps({'status': 'downloading', 'started_utc': started,
                                  'destination': str(DEST), 'expected_bytes': before['Size']}, indent=2) + '\n')
t0 = time.perf_counter()
subprocess.run(command, check=True)
hashes = {'sha256': hashlib.sha256(), 'md5': hashlib.md5(), 'sha1': hashlib.sha1()}
with DEST.open('rb') as source:
    for block in iter(lambda: source.read(8 << 20), b''):
        for h in hashes.values():
            h.update(block)
actual = {name: h.hexdigest() for name, h in hashes.items()}
after = metadata()
assert before['Hashes'] == after['Hashes'] and before['ModTime'] == after['ModTime']
assert actual == before['Hashes'] and DEST.stat().st_size == before['Size']
receipt = {'status': 'complete', 'started_utc': started,
           'completed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
           'seconds_including_verification': time.perf_counter() - t0,
           'drive_url': 'https://drive.google.com/file/d/' + FILE_ID + '/view',
           'drive_file_id': FILE_ID, 'remote': REMOTE, 'local_file': str(DEST),
           'bytes': DEST.stat().st_size, 'hashes': actual,
           'drive_identity_unchanged': True, 'provider_hashes_match': True,
           'command': command, 'download_log': str(WORK / 'download.log')}
status_path.write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(receipt, indent=2), flush=True)
