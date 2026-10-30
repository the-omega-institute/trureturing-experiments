from itertools import product
from collections import Counter
import time,json,sys

def pos_le(a,i,j):
    if i==j:return True
    if i>j:return False
    return all(a[k]>a[j] for k in range(i,j))
def le(a,b):
    if any(x>y for x,y in zip(a,b)): return False
    n=len(a)
    for j in range(n):
      for i in range(j):
        if pos_le(a,i,j) and a[i]-a[j]>b[i]-b[j]: return False
    return True

def parent(a):
 n=len(a); p=[-1]*n
 for i in range(n):
  for j in range(i+1,n):
   if a[j]<a[i]: p[i]=j; break
 return p

def ht(a,b):
 p=parent(a); s=0
 for i,(x,y) in enumerate(zip(a,b)):
  d=y-x
  if p[i]>=0:d-=b[p[i]]-a[p[i]]
  if d<0: raise AssertionError((a,b,p,d))
  s+=d
 return s

def poly(b):
 c=Counter(); n=len(b); cnt=0
 for a in product(*(range(x+1) for x in b)):
  if le(a,b): c[ht(a,b)]+=1;cnt+=1
 return [c[i] for i in range(max(c)+1)] if c else []

def search(n,maxv):
 t=time.monotonic(); checks=0; polys={};
 for b in product(range(maxv+1),repeat=n):
  pb=poly(b); rb=b[::-1]; pr=poly(rb); checks+=1
  if pb!=pr:
   return {'status':'counterexample','n':n,'b':b,'rev':rb,'f_b':pb,'f_rev':pr,'checked':checks,'elapsed':time.monotonic()-t}
 return {'status':'no_counterexample','n':n,'maxv':maxv,'checked':checks,'elapsed':time.monotonic()-t}
if __name__=='__main__':
 for n,m in [(1,12),(2,10),(3,7),(4,5),(5,3)]:
  r=search(n,m);print(json.dumps(r),flush=True)
  if r['status']=='counterexample':break
