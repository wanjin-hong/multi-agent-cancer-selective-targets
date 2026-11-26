"""Download GDC STAR-count files one by one (resilient to /data tarball
stalls). Skips file_ids already extracted under gdc_subset/.
Each file saved as gdc_subset/<file_id>/<file_name>.
"""
import json, os, time, requests

DL = '/home/hatch/workspace/research_notes/multi-cancer-selective-targets-20261120/neuroblastoma/downloads'
subset = json.load(open(os.path.join(DL, 'gdc_subset_files.json')))
all_hits = subset['nbl'] + subset['pooled']
done = set(os.listdir(os.path.join(DL, 'gdc_subset'))) if os.path.isdir(os.path.join(DL, 'gdc_subset')) else set()
print('already have:', len(done), 'of', len(all_hits), flush=True)
n_ok = n_skip = n_fail = 0
fails = []
for h in all_hits:
    fid = h['file_id']
    if fid in done:
        n_skip += 1
        continue
    fname = h['file_name']
    outdir = os.path.join(DL, 'gdc_subset', fid)
    out = os.path.join(outdir, fname)
    ok = False
    for attempt in range(4):
        try:
            r = requests.get(f'https://api.gdc.cancer.gov/data/{fid}', timeout=300, stream=True)
            if r.status_code == 200:
                os.makedirs(outdir, exist_ok=True)
                with open(out, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        f.write(chunk)
                if os.path.getsize(out) > 100000:
                    ok = True
                    break
            else:
                print('HTTP', r.status_code, 'for', fid, flush=True)
        except Exception as e:
            print('attempt', attempt, 'failed for', fid, ':', e, flush=True)
        time.sleep(5)
    if ok:
        n_ok += 1
    else:
        n_fail += 1
        fails.append(fid)
    if (n_ok + n_fail) % 20 == 0:
        print(f'progress: ok {n_ok} fail {n_fail} skip {n_skip}', flush=True)
print('DONE ok', n_ok, 'fail', n_fail, 'skip', n_skip, flush=True)
print('failed ids:', fails, flush=True)
json.dump(fails, open(os.path.join(DL, 'gdc_failed_ids.json'), 'w'))
