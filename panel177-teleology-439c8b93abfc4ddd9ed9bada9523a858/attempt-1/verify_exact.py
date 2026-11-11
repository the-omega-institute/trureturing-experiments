"""Solver-independent exact rational verifier of the retained local certificate.
No floating-point computation or package beyond the Python standard library.
"""
from fractions import Fraction as F
from pathlib import Path
import itertools,json
ROOT=Path(__file__).resolve().parent
d=json.loads((ROOT/'exact_active8_simplex.json').read_text());L=d['L'];words=list(itertools.product((0,1),repeat=L));first=lambda w:next((i+1 for i,x in enumerate(w) if not x),0)
reward={w:[F(0)]*(L+1) for w in words};beta=[F(0)]*L;K=C0=W=F(0)
for tag,word,v in d['cut']['tokens']:
 w=tuple(map(int,word));n=first(w);v=F(v);assert v>=0;K+=v
 if tag in ('A','P','T'):reward[w][0]+=v;C0+=v
 if tag in ('C','P'):reward[w][n]+=v
 if tag=='A':beta[n-1]-=v
 if tag=='C':beta[n-1]+=v
 anchor=1 if w[0]==0 else 2 if w[1]==0 else 0
 if (tag=='A' and anchor==0) or (tag=='C' and anchor==n) or tag in ('P','T'):W+=v
B=F(0);full_schedules=[()]+[(i,) for i in range(L)]+list(itertools.combinations(range(L),2));cellcount=0
for S,record in zip(full_schedules,d['cut']['complete_N_cell_maxima']):
 assert record['schedule']==[i+1 for i in S];cells={}
 for w in words:
  cell=tuple(w[i] for i in S);r=cells.setdefault(cell,[F(0)]*(L+1))
  for z,v in enumerate(reward[w]):r[z]+=v
 maxima=[max(cells[c]) for c in sorted(cells)];assert maxima==list(map(F,record['maxima']));B=max(B,sum(maxima));cellcount+=len(cells)
assert len(d['cut']['complete_N_cell_maxima'])==len(full_schedules)
assert B==W==54 and K==138 and C0==51 and beta==list(map(F,d['cut']['beta']))
base=[F(x,64) for x in d['base_b64']];step=F(d['simplex_step']);directions=d['directions'];assert len(d['vertices'])==8
checks=0;tablecount=0
for vi,vertex in enumerate(d['vertices']):
 b=list(map(F,vertex['b']));e=F(vertex['e']);ids=[] if vi==0 else directions[vi-1]
 assert b==[x+step*int(i+1 in ids) for i,x in enumerate(base)]
 assert e==(sum(x*y for x,y in zip(beta,b))-3)/K
 assert b[0]<=1 and all(b[i]>=b[i+1] for i in range(7)) and b[1]==b[2] and b[-1]>e>0
 assert 9*e<=sum(b[2:]) and b[0]+sum(b[2:])-9*e<=1
 tables=vertex['tables'];mass=F(0)
 for t in tables:
  S=t['schedule'];assert len(S)<=2 and S==sorted(set(S)) and all(1<=i<=L for i in S)
  weight=F(t['mass']);assert weight>=0;mass+=weight
  assert len(t['cells'])==(1 if not S else 2 if len(S)==1 else 3)
  for c in t['cells']:
   assert all(0<=int(z)<=L and F(v)>=0 for z,v in c.items());assert sum(map(F,c.values()))==weight
 assert mass==1;tablecount+=len(tables)
 minima=[None]*L;wrong=[F(0)]*L
 for w in words:
  law=[F(0)]*(L+1)
  for t in tables:
   S=t['schedule'];c=0 if not S else w[S[0]-1] if len(S)==1 else 0 if w[S[0]-1]==0 else 1+w[S[1]-1]
   for z,v in t['cells'][c].items():law[int(z)]+=F(v)
  assert sum(law)==1;n=first(w)
  if n:
   assert law[0]>=1-b[n-1]-e;assert law[n]>=b[n-1]-e;assert law[0]+law[n]>=1-e;checks+=3
   minima[n-1]=law[n] if minima[n-1] is None else min(minima[n-1],law[n]);wrong[n-1]=max(wrong[n-1],1-law[0]-law[n])
  else:
   assert law[0]>=1-e;checks+=1
   assert {str(z):str(v) for z,v in enumerate(law) if v}==vertex['positive_law']
 assert minima==list(map(F,vertex['correct_minima'])) and wrong==list(map(F,vertex['wrong_maxima']))
assert checks==6128
result={'status':'exact_pass','vertices':8,'free_amplitude_dimension':7,'word_laws':2048,'event_inequalities':checks,'positive_mass_schedule_rows':tablecount,'full_N_schedules':len(full_schedules),'full_N_cells':cellcount,'N_support':str(B),'N_anchor':str(W),'score_K':str(K),'scope':'ordinary rational finite certificate only; no global selector, strict separator or Lean claim'}
(ROOT/'exact_verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
