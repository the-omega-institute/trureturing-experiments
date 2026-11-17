import json,time
from pathlib import Path
import sympy as s
P=Path(__file__).parent;b=s.symbols('b');t=time.perf_counter();out=[]
for i in range(9):
 d=json.loads((P/('support_factor_'+str(i)+'.json')).read_text());p=s.Poly.from_dict({(int(k[1]),):s.Rational(c) for k,c in d['terms']},b,domain=s.QQ)
 ss=time.perf_counter();ints=p.intervals(eps=s.Rational(1,10**12));out.append({'index':i,'degree':p.degree(),'real_roots':len(ints),'seconds':time.perf_counter()-ss,'intervals':[[str(a),str(z),m] for (a,z),m in ints]})
 (P/'support_isolations.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out[-1].items() if k!='intervals'}),flush=True)
print(json.dumps({'total_seconds':time.perf_counter()-t,'total_real_roots':sum(x['real_roots'] for x in out)}))
