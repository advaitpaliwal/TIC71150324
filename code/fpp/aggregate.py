from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
D=Path('work/deep/71150324/fpp')
S=Path('work/hunt7/fpp/71150324')
original=json.loads((S/'out/res_71150324_0.json').read_text())
rows=[]; scenarios=[]
for k in range(10):
 r=json.loads((D/f'out/res_71150324_{k}.json').read_text())
 for field in ['P','T0','depth_ppm','err_ppm','nbins']:
  assert r[field]==original[field], (k,field,r[field],original[field])
 assert r['N']==1000000 and r['seed']==1000+k and len(r['dropped'])==7
 assert r['stars_before']==45 and r['stars_after']==38 and r['pixel_shapes']==[[38,2]]*3
 p=pd.read_csv(D/f'out/probs_71150324_{k}.csv',dtype={'ID':str})
 assert not p.ID.isin(r['dropped']).any()
 assert np.isclose(p.prob.sum(),1)
 # Installed TRICERATOPS uses row positions 0,3,9 = TP,PTP,DTP, not STP.
 assert np.isclose(1-p.iloc[[0,3,9]].prob.sum(),r['FPP'])
 assert np.isclose(p.loc[p.scenario.str.startswith('N'),'prob'].sum(),r['NFPP'])
 row={key:r[key] for key in ['run','seed','N','depth_ppm','err_ppm','FPP','NFPP']}
 for name in ['TP','PTP','STP','NEB','NEBx2P','NTP']:
  row[name]=float(p.loc[p.scenario==name,'prob'].sum())
 row['EBx2P_sum']=float(p.loc[p.scenario.str.endswith('EBx2P'),'prob'].sum())
 row['any_NEB']=row['NEB']+row['NEBx2P']
 rows.append(row)
 for _,s in p.iterrows(): scenarios.append({'run':k,'ID':s.ID,'scenario':s.scenario,'prob':s.prob})
a=pd.DataFrame(rows);a.to_csv(D/'runs.csv',index=False)
b=pd.DataFrame(scenarios); b.groupby(['ID','scenario']).prob.agg(['mean','std','min','max']).reset_index().to_csv(D/'scenario_summary.csv',index=False)
summary={'runs':10,'N_per_run':1000000,'sample_sd_ddof':1,'removed_stars':7,'remaining_stars':38,'depth_ppm':original['depth_ppm'],'err_ppm':original['err_ppm'],'statistics':{}}
for name in ['FPP','NFPP','TP','PTP','STP','EBx2P_sum','NEB','NEBx2P','NTP','any_NEB']:
 v=a[name]; summary['statistics'][name]={'mean':float(v.mean()),'sample_sd':float(v.std(ddof=1)),'min':float(v.min()),'max':float(v.max())}
summary['input_sha256']={}
for f in ['71150324_gaia_background.csv','out/target_71150324.pkl','molusc_out/molusc_for_trc.csv']:
 h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
 assert h(S/f)==h(D/f)
 summary['input_sha256'][f]=h(D/f)
(D/'joint_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
print(b[b.scenario.str.startswith('N')].groupby(['ID','scenario']).prob.agg(['mean','std','min','max']).to_string())
