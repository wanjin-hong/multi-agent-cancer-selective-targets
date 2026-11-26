"""Retry the 56 failed GDC downloads with curl (urllib kept hitting
IncompleteRead on these). Verifies each file against its expected
file_size from the GDC query metadata; retries up to 5 times.
"""
import json, os, subprocess

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/neuroblastoma/downloads'
subset = json.load(open(os.path.join(DL, 'gdc_subset_files.json')))
all_hits = {h['file_id']: h for h in subset['nbl'] + subset['pooled']}
fails = json.load(open(os.path.join(DL, 'gdc_failed_ids_threaded.json')))
print('to retry:', len(fails))
still = []
for fid in fails:
    h = all_hits[fid]
    fname = h['file_name']
    expected = h['file_size']
    outdir = os.path.join(DL, 'gdc_subset', fid)
    out = os.path.join(outdir, fname)
    os.makedirs(outdir, exist_ok=True)
    ok = False
    for attempt in range(5):
        r = subprocess.run(['curl', '-s', '--retry', '2', '--max-time', '300',
                            f'https://api.gdc.cancer.gov/data/{fid}', '-o', out],
                           capture_output=True, text=True)
        if os.path.exists(out) and os.path.getsize(out) >= expected - 1000:
            ok = True
            break
    print(fid, 'OK' if ok else 'FAIL', os.path.getsize(out) if os.path.exists(out) else 0, 'expected', expected, flush=True)
    if not ok:
        still.append(fid)
print('DONE ok', len(fails) - len(still), 'still failing', len(still))
print('still failing:', still)
json.dump(still, open(os.path.join(DL, 'gdc_failed_ids_curl.json'), 'w'))
