# Fill rnaas_rev1.md placeholders from summary_full.json / summary_sel.json / perstar_nfpp_full.csv.
import json, csv
O='~/Scratch/planet-hunt/revision-82162/'
F=json.load(open(O+'fpp/summary_full.json')); S=json.load(open(O+'fpp/summary_sel.json'))
sh={r['tic']:float(r['share_mean']) for r in csv.DictReader(open(O+'fpp/perstar_nfpp_full.csv'))}
def sci(x):
    m,e=f'{x:.1e}'.split('e'); return f'{m} × 10⁻{int(e[1:])}'.replace('⁻','⁻'+'') if False else f'{m} × 10^{int(e)}'
sup=str.maketrans('-0123456789','⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
def sci(x): m,e=f'{x:.1e}'.split('e'); return f'{m} × 10'+str(int(e)).translate(sup)
v={'@FPP_FULL@':f"{F['FPP_mean']:.3f} ± {F['FPP_sd']:.3f}", '@NFPP_FULL@':f"{F['NFPP_mean']:.4f}",
   '@SHARE_20@':f"{100*sh['71150320']:.0f}%", '@NFPP_SEL@':sci(S['NFPP_mean']), '@FPP_SEL@':f"{S['FPP_mean']:.3f} ± {S['FPP_sd']:.3f}",
   '@STP_FULL@':f"{F['scenarios']['STP']:.3f}"}
p=O+'final/rnaas_rev1.md'; s=open(p).read()
for k,x in v.items(): s=s.replace(k,x)
assert '@' not in s.replace('@gmail',''); open(p,'w').write(s); print(v)
