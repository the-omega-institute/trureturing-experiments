"""Exact certificate of one local active-eight class-equality simplex.
Success (preregistered before runs): all eight rational stopped lifts obey all
cell equalities and all 766 event demands, and one complete N score gives
an exactly tight common lower plane. Fail at any nonzero rational residual.
These finite identities do not imply a global shared score selector.
"""
from pathlib import Path
from fractions import Fraction as F
from math import lcm,gcd
import itertools,json
import numpy as np
from search_complete import model,solve
ROOT=Path(__file__).resolve().parent
base=json.loads((ROOT/'frontier-candidate-0.json').read_text());B=[F(x,64) for x in base['a64']];L=8
score=[(tag,tuple(map(int,w)),F(v)*138) for tag,w,v in base['N_y']]
assert all(v.denominator==1 for _,_,v in score)
words=list(itertools.product((0,1),repeat=L));first=lambda w:next((i+1 for i,x in enumerate(w) if not x),0)
f={w:[F(0)]*9 for w in words};beta=[F(0)]*8;K=F(0);C0=F(0);W=F(0)
for tag,w,v in score:
 n=first(w); K+=v
 if tag in ('A','P','T'):f[w][0]+=v;C0+=v
 if tag in ('C','P'):f[w][n]+=v
 if tag=='C':beta[n-1]+=v
 if tag=='A':beta[n-1]-=v
 anchor=first(w[:2]) or 0
 if (tag=='A' and anchor==0) or (tag=='C' and anchor==n) or tag in ('P','T'):W+=v
cellmax=[]
for S in [()]+[(i,) for i in range(L)]+list(itertools.combinations(range(L),2)):
 cells={}
 for w in words:
  c=tuple(w[i] for i in S)
  if c not in cells:cells[c]=[F(0)]*9
  cells[c]=[x+y for x,y in zip(cells[c],f[w])]
 vals=[max(cells[c]) for c in sorted(cells)]
 cellmax.append({'schedule':[i+1 for i in S],'maxima':[str(x) for x in vals]})
BN=max(sum(F(x) for x in s['maxima']) for s in cellmax)
assert beta==[F(x) for x in [3,0,12,15,15,15,15,12]]
assert K==138 and BN==W
# One simplex direction for each free amplitude on b2=b3.
dirs=[[i] for i in [0,3,4,5,6,7]];dirs.insert(1,[1,2])
m=model(8,'O');vertices=[]
for vi,ids in enumerate([[]]+dirs):
 b=[x+F(1,4096)*int(i in ids) for i,x in enumerate(B)]
 e=(sum(x*y for x,y in zip(beta,b))+C0-BN)/K
 assert b[-1]>e>0 and all(b[i]>=b[i+1] for i in range(7)) and b[1]==b[2]
 assert 9*e<=sum(b[2:]) and b[0]+sum(b[2:])-9*e<=1
 r=solve(m,np.array([float(x) for x in b]));assert abs(r.fun-float(e))<1e-8
 support=np.flatnonzero(r.x>1e-9)
 eqmat=m['E'][:,support].toarray();eqrhs=list(map(F,m['eqrhs']))
 d=[]
 for n in m['ns']:d.extend([1-b[n-1],b[n-1],F(1)] if n else [F(1)])
 active=np.flatnonzero(r.ineqlin.residual<1e-8)
 ubmat=m['A'][active,:][:,support].toarray();ubrhs=[-d[i] for i in active]
 mat=[];rhs=[]
 for row,v in list(zip(eqmat,eqrhs))+list(zip(ubmat,ubrhs)):
  if np.any(row):mat.append([int(v) for v in row]);rhs.append(v)
 basis={}
 for row,v in zip(mat,rhs):
  work={i:F(z) for i,z in enumerate(row) if z};value=v
  while work:
   j=min(work)
   if j not in basis:
    z=work[j];work={k:x/z for k,x in work.items()};value/=z;basis[j]=(work,value);break
   old,oldvalue=basis[j];z=work[j]
   for k,xold in old.items():
    work[k]=work.get(k,F(0))-z*xold
    if not work[k]:work.pop(k)
   value-=z*oldvalue
  else:assert value==0
 assert len(basis)==len(support)
 sol=[F(0)]*len(support)
 for j in sorted(basis,reverse=True):
  row,v=basis[j];sol[j]=v-sum(x*sol[k] for k,x in row.items() if k!=j)
 x=[F(0)]*len(r.x)
 for col,v in zip(support,sol):x[int(col)]=v
 assert x[-1]==e
 tables=[];mass=F(0)
 for s in m['catalog']:
  weight=x[s['mass_col']];cells=[[x[c] for c in cols] for cols in s['cell_cols']]
  assert all(v>=0 for c in cells for v in c)
  assert all(sum(c)==weight for c in cells)
  mass+=weight
  if weight:tables.append({'schedule':s['schedule'],'mass':str(weight),'cells':[{str(z):str(v) for z,v in enumerate(c) if v} for c in cells]})
 assert mass==1
 minima=[None]*8;maxwrong=[F(0)]*8;rows=0;min_slack=None;positive=None
 for w in words:
  law=[F(0)]*9
  for t in tables:
   S=[i-1 for i in t['schedule']]
   ci=0 if not S else w[S[0]] if len(S)==1 else 0 if not w[S[0]] else 1+w[S[1]]
   for z,v in t['cells'][ci].items():law[int(z)]+=F(v)
  assert sum(law)==1;n=first(w)
  if n:
   slack=[law[0]-(1-b[n-1]-e),law[n]-(b[n-1]-e),law[0]+law[n]-(1-e)];rows+=3
   assert min(slack)>=0
   min_slack=min(slack) if min_slack is None else min(min_slack,*slack)
   minima[n-1]=law[n] if minima[n-1] is None else min(minima[n-1],law[n]);maxwrong[n-1]=max(maxwrong[n-1],1-law[0]-law[n])
  else:
   assert law[0]>=1-e;positive={str(z):str(v) for z,v in enumerate(law) if v};rows+=1
 assert rows==766
 vertices.append({'index':vi,'b':[str(t) for t in b],'e':str(e),'tables':tables,'correct_minima':[str(v) for v in minima],'wrong_maxima':[str(v) for v in maxwrong],'positive_law':positive,'event_rows_checked':rows,'min_slack':str(min_slack)})
# Integer form for vertex0's complete table.
tables=vertices[0]['tables'];den=lcm(*(F(t['mass']).denominator for t in tables),*(F(v).denominator for t in tables for c in t['cells'] for v in c.values()))
intables=[{'schedule':t['schedule'],'mass':int(F(t['mass'])*den),'cells':[{z:int(F(v)*den) for z,v in c.items()} for c in t['cells']]} for t in tables]
output={'schema':'panel177-teleology-exact-active8-local-simplex-v1','L':8,'base_b64':base['a64'],'simplex_step':'1/4096','directions':[[i+1 for i in ids] for ids in dirs],'cut':{'tokens':[[tag,''.join(map(str,w)),int(v)] for tag,w,v in score],'beta':[str(x) for x in beta],'K':str(K),'C0':str(C0),'anchor':str(W),'BN':str(BN),'complete_N_cell_maxima':cellmax},'vertices':vertices,'primary_denominator':den,'primary_integer_tables':intables,'proof_status':'exact rational arithmetic, no Lean/kernel claim; local simplex only','total_event_rows_checked':8*766,'primary_cut_primitive':'b1+4b3+5(b4+b5+b6+b7)+4b8 <= 1+46e'}
(ROOT/'exact_active8_simplex.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps({'verified_vertices':8,'event_rows':8*766,'tables_per_vertex':[len(v['tables']) for v in vertices],'BN':str(BN),'W':str(W),'C0':str(C0),'K':str(K),'base_error':vertices[0]['e'],'denominator':den,'beta':output['cut']['beta']}))
