"""Zhang (arXiv:2303.01024v1) closing question: are two T2-threshold k-hypergraphs with the same independence polynomial
isomorphic? Family (k = 3, N >= 5): H_N = all triples of [N] except 123, 124, 125; K_N = all except 123, 124, 134.
Checks for N = 5..8: both satisfy (T2) with explicit labels (H: c = 0,0,1,1,1,2,...; tau = 1. K: d = 0,1,1,1,3,...;
tau = 2), equal independence polynomials 1 + N x + C(N,2) x^2 + 3 x^3, and no vertex bijection maps H_N onto K_N.
"""
from itertools import combinations, permutations
def hyper(N, missing):
    V=range(1,N+1); return {frozenset(t) for t in combinations(V,3)} - {frozenset(m) for m in missing}
def t2(N,E,c,tau):
    return all((frozenset(t) in E) == (sum(c[v] for v in t) > tau) for t in combinations(range(1,N+1),3))
def indpoly(N,E):
    co=[0]*(N+1)
    for r in range(N+1):
        for S in combinations(range(1,N+1),r):
            s=set(S)
            if not any(e<=s for e in E): co[r]+=1
    return co
def iso(N,E1,E2):
    if len(E1)!=len(E2): return False
    for p in permutations(range(1,N+1)):
        m={i+1:p[i] for i in range(N)}
        if {frozenset(m[v] for v in e) for e in E1}==E2: return True
    return False
ok=True
for N in range(5,9):
    H=hyper(N,[(1,2,3),(1,2,4),(1,2,5)]); K=hyper(N,[(1,2,3),(1,2,4),(1,3,4)])
    cH={1:0,2:0,3:1,4:1,5:1,**{i:2 for i in range(6,N+1)}}; cK={1:0,2:1,3:1,4:1,**{i:3 for i in range(5,N+1)}}
    pH=indpoly(N,H); pK=indpoly(N,K)
    a,b,c2,d=t2(N,H,cH,1),t2(N,K,cK,2),pH==pK,iso(N,H,K); ok = ok and a and b and c2 and not d
    print(N, 'T2(H)=',a,'T2(K)=',b,'poly equal=',c2,pH,'isomorphic=',d)
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
