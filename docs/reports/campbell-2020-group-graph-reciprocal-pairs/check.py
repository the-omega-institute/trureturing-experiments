"""Exact check of a counterexample family to Campbell's conjecture (arXiv:2005.13566) on irreducible
group-graph reciprocal pairs: Gamma = 3K_a + 2K_b (disjoint), G = (S_a^3 x S_b^2) x| S_3 where S_3 permutes the
three K_a naturally and the two K_b through the sign. Checks P_{Gamma,G}(x) = (-1)^n F_G(-x) by enumerating G."""
import sys
from itertools import permutations
from functools import lru_cache
def poly_add(p,q):
    r=[0]*max(len(p),len(q))
    for i,c in enumerate(p): r[i]+=c
    for i,c in enumerate(q): r[i]+=c
    return r
def poly_sub(p,q): return poly_add(p,[-c for c in q])
@lru_cache(maxsize=None)
def chrom(nv,edges):
    # chromatic polynomial (coefficient list, index = degree) of graph on vertices 0..nv-1 with frozenset edges
    if not edges:
        r=[0]*(nv+1); r[nv]=1; return tuple(r)
    e=min(edges); u,v=e
    rest=frozenset(edges-{e})
    a=chrom(nv,rest)
    # contract v into u, relabel
    def m(x): 
        x=u if x==v else x
        return x-1 if x>v else x
    ce=frozenset(tuple(sorted((m(x),m(y)))) for x,y in rest if m(x)!=m(y))
    b=chrom(nv-1,ce)
    return tuple(poly_sub(list(a),list(b)))
def check(a,b):
    blocks=[list(range(i*a,(i+1)*a)) for i in range(3)]+[list(range(3*a+j*b,3*a+(j+1)*b)) for j in range(2)]
    n=3*a+2*b
    E=set()
    for bl in blocks:
        for i in bl:
            for j in bl:
                if i<j: E.add((i,j))
    # build all elements: choose internal perms in each block and a top perm h in S3
    def sign(h): 
        s=1
        for i in range(3):
            for j in range(i+1,3):
                if h[i]>h[j]: s=-s
        return s
    F={}; P={}
    from collections import Counter
    Fc=Counter(); Pc=Counter()
    pa=list(permutations(range(a))); pb=list(permutations(range(b)))
    for h in permutations(range(3)):
        hb=(0,1) if sign(h)==1 else (1,0)
        for t1 in pa:
          for t2 in pa:
            for t3 in pa:
              ts=(t1,t2,t3)
              for u1 in pb:
                for u2 in pb:
                  us=(u1,u2)
                  g=[0]*n
                  for i in range(3):
                      for k in range(a): g[blocks[i][k]]=blocks[h[i]][ts[i][k]]
                  for j in range(2):
                      for k in range(b): g[blocks[3+j][k]]=blocks[3+hb[j]][us[j][k]]
                  # cycles
                  seen=[False]*n; cyc=[]
                  for s in range(n):
                      if not seen[s]:
                          c=[]; x=s
                          while not seen[x]: seen[x]=True; c.append(x); x=g[x]
                          cyc.append(c)
                  Fc[len(cyc)]+=1
                  idx={}
                  for ci,c in enumerate(cyc):
                      for x in c: idx[x]=ci
                  loop=any(idx[x]==idx[y] for x,y in E)
                  if loop: continue
                  qe=frozenset(tuple(sorted((idx[x],idx[y]))) for x,y in E)
                  Pc[(len(cyc),qe)]+=1
    Fp=[0]*(n+1)
    for k,c in Fc.items(): Fp[k]+=c
    Pp=[0]*(n+1)
    for (nv,qe),c in Pc.items():
        cp=chrom(nv,qe)
        for i,co in enumerate(cp): Pp[i]+=c*co
    order=sum(Fc.values())
    rhs=[((-1)**n)*((-1)**i)*Fp[i] for i in range(n+1)]   # (-1)^n F(-x)
    ok=(Pp==rhs)
    print(f"a={a} b={b} n={n} |G|={order} reciprocal={ok}")
    print("  F_G coeffs (deg 0..n):",Fp)
    return ok
ok=all(check(a,b) for a,b in [(2,4),(4,2)])
print('ALL_OK' if ok else 'FAIL'); sys.exit(0 if ok else 1)
