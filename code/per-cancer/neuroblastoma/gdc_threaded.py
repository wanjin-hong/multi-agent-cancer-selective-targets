"""Threaded GDC STAR-count download (6 workers). Skips file_ids already
extracted under gdc_subset/. Same layout as gdc_perfile.py.
"""
import json, os, time, requests
from concurrent.futures import ThreadPoolExecutor

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/neuroblastoma/downloads'
LOG = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/neuroblastoma/gdc_download_threaded.log'
subset = json.load(open(os.path.join(DL, 'gdc_subset_files.json')))
all_hits = subset['nbl'] + subset['pooled']
done = set(os.listdir(os.path.join(DL, 'gdc_subset'))) if os.path.isdir(os.path.join(DL, 'gdc_subset')) else set()
todo = [h for h in all_hits if h['file_id'] not in done]
print(f'already have: {len(done)} of {len(all_hits)}; to download: {len(todo)}', flush=True)

def fetch(h):
    fid = h['file_id']
    fname = h['file_name']
    outdir = os.path.join(DL, 'gdc_subset', fid)
    out = os.path.join(outdir, fname)
    for attempt in range(4):
        try:
            r = requests.get(f'https://api.gdc.cancer.gov/data/{fid}', timeout=300, stream=True)
            if r.status_code == 200:
                os.makedirs(outdir, exist_ok=True)
                with open(out, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        f.write(chunk)
                if os.path.getsize(out) > 100000:
                    return (fid, True, '')
            else:
                last = f'HTTP {r.status_code}'
        except Exception as e:
            last = str(e)[:120]
        time.sleep(5)
    return (fid, False, last)

n_ok = n_fail = 0
fails = []
with open(LOG, 'a') as log:
    with ThreadPoolExecutor(max_workers=6) as ex:
        for i, (fid, ok, err) in enumerate(ex.map(fetch, todo)):
            if ok:
                n_ok += 1
            else:
                n_fail += 1
                fails.append(fid)
                log.write(f'FAIL {fid} {err}\n')
                log.flush()
            if (i + 1) % 20 == 0:
                print(f'progress: ok {n_ok} fail {n_fail} of {len(todo)}', flush=True)
print('DONE ok', n_ok, 'fail', n_fail, 'skip', len(done), flush=True)
print('failed ids:', fails, flush=True)
json.dump(fails, open(os.path.join(DL, 'gdc_failed_ids_threaded.json'), 'w'))
