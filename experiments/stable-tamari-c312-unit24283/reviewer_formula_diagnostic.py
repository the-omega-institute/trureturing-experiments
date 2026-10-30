from fractions import Fraction as F
import json
from pathlib import Path

def mul(p,q):
 r={}
 for u,a in p.items():
  for v,b in q.items():
   w=tuple(x+y for x,y in zip(u,v));r[w]=r.get(w,0)+a*b
 return {w:a for w,a in r.items() if a}
def add(*ps):
 r={}
 for p in ps:
  for w,a in p.items():r[w]=r.get(w,0)+a
 return {w:a for w,a in r.items() if a}
one={(0,0,0):1};x={(1,0,0):1};y={(0,1,0):1}
def neg(p):return {w:-a for w,a in p.items()}
def sub(p,q):return add(p,neg(q))
txy={(1,1,1):1};tx={(1,0,1):1};ty={(0,1,1):1}
# Common denominator: (1-xy)(1-x)(1-y)(1-tx)(1-ty)(1-txy).
N=sub(add(mul(sub(one,x),sub(one,txy)),mul(mul(x,sub(one,y)),sub(one,ty))),add(mul(sub(one,y),sub(one,txy)),mul(mul(y,sub(one,x)),sub(one,tx))))
assert N=={}
xv,yv,tv=F(1,5),F(1,7),F(2)
def W0(x,y,t):return 1/((1-x*y)*(1-y)*(1-t*x)*(1-t*y))
def W1(x,y,t):return x/((1-x*y)*(1-x)*(1-t*x)*(1-t*x*y))
assert W0(xv,yv,tv)+W1(xv,yv,tv)==W0(yv,xv,tv)+W1(yv,xv,tv)
assert W1(xv,yv,tv) not in (W0(yv,xv,tv), W1(yv,xv,tv))
r={'scope':'Theorem 4.17 n=2 diagnostic only; no conjecture search','common_denominator_difference':N,'reflected_intervals_h_11':[[2,2],[1,2]],'endpoint_collision':True,'term_permutation_fails_at':{'x':'1/5','y':'1/7','t':'2'},'W1':str(W1(xv,yv,tv)),'reflected_W0':str(W0(yv,xv,tv)),'reflected_W1':str(W1(yv,xv,tv)),'full_sum_identity_checked':'exact integer polynomial expansion'}
with (Path(__file__).resolve().parent / 'reviewer_formula_diagnostic.json').open('w') as f:json.dump(r,f,indent=2)
print(json.dumps(r))
