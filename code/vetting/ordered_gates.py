# Glue over hunt5a/raw_gate.py summaries; never invoke later gates after a failure.
from pathlib import Path
import json, math
r=Path(__file__).resolve().parent
rows=[]
for f in sorted(r.glob('raw_gate_*.json')):
 x=json.loads(f.read_text()); failures=[]; per=x['per_sector']; combined=x['combined']
 def z(a,b): return abs(a['depth_ppm']-b['depth_ppm'])/math.hypot(a['err_ppm'],b['err_ppm'])
 if not per or not combined: failures.append('No measurable raw/corrected events')
 for key in x.get('available',[]):
  if key not in per:failures.append(key+': product has data but no measurable event depth')
 comparisons={}
 for key,v in per.items():
  if not v.get('raw') or not v.get('corrected'):
   failures.append(key+': missing finite depth');continue
  zz=z(v['raw'],v['corrected']); comparisons[key]=zz
  if not math.isfinite(zz) or zz>3:failures.append(f'{key}: raw/corrected disagree {zz:.4f} sigma')
 raw_snr={}
 for prod,v in combined.items():
  if not v.get('raw'): failures.append(prod+': no finite raw depth');continue
  raw=v['raw'];raw_snr[prod]=raw['depth_ppm']/raw['err_ppm']
  if not raw['depth_ppm']>0:failures.append(prod+': combined raw nonpositive')
 if not any(v>=3 for v in raw_snr.values()):failures.append('No product combined raw >=3 sigma')
 row={'tic':x['tic'],'P':x['P'],'a':{'status':'fail' if failures else 'pass','failures':failures,'raw_snr':raw_snr,'raw_corrected_sigma':comparisons},'b':{'status':'not run'},'c':{'status':'not run'},'d':{'status':'not run'}}
 if not failures:
  failures=[];cross={}
  sectors=sorted({int(k.split('_')[1]) for k in per})
  paired=[s for s in sectors if f'spoc_{s}' in per and f'qlp_{s}' in per]
  qualified=[s for s in paired if all(per[f'{p}_{s}'][label]['events']>=3 for p in ['spoc','qlp'] for label in ['raw','corrected'])]
  if len(qualified)<2:failures.append('Fewer than two paired sectors with >=3 events in both products')
  values=[per[f'{p}_{s}']['corrected'] for s in paired for p in ['spoc','qlp']]
  weight=sum(1/v['err_ppm']**2 for v in values)
  reference={'depth_ppm':sum(v['depth_ppm']/v['err_ppm']**2 for v in values)/weight,'err_ppm':1/math.sqrt(weight)} if weight else None
  unpaired={}
  for s in sectors:
   a=per.get(f'spoc_{s}');b=per.get(f'qlp_{s}')
   if not a or not b:
    v=(a or b)['corrected'];zz=z(v,reference) if reference else float('nan');unpaired[str(s)]={'depth_ppm':v['depth_ppm'],'reference_sigma':zz}
    if v['depth_ppm']<=0 or not math.isfinite(zz) or zz>3:failures.append(f'sector {s}: unpaired nonpositive or inconsistent {zz:.4f} sigma')
    continue
   zz=z(a['corrected'],b['corrected']);cross[str(s)]=zz
   if not math.isfinite(zz) or zz>3:failures.append(f'sector {s}: products disagree {zz:.4f} sigma')
  row['b']={'status':'fail' if failures else 'pass','failures':failures,'corrected_cross_sigma':cross,'unpaired_coverage':unpaired,'paired_reference':reference,'events':{k:v['corrected']['events'] for k,v in per.items()}}
 if row['b']['status']=='pass':
  loc=r/f"localize/{x['tic']}/gate_c.json"
  if loc.exists():row['c']=json.loads(loc.read_text())
 rows.append(row)
(r/'ordered_gates.json').write_text(json.dumps(rows,indent=2))
for row in rows: print(row['tic'],row['P'],[(k,row[k]['status']) for k in 'abcd'],row['a']['failures'],row['b'].get('failures',[]))
