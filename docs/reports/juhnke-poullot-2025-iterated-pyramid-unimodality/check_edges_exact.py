import random, itertools
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog
from sympy import Matrix
random.seed(3)
R=[tuple(Fr(random.randint(-20,20),random.randint(1,5)) for _ in range(4)) for _ in range(9)]
def hyperplane(pts):
    # normal n with n.(p_i - p_0) = 0 for i=1..3 (exact nullspace)
    M=Matrix([[pts[i][t]-pts[0][t] for t in range(4)] for i in range(1,4)])
    ns=M.nullspace()
    if len(ns)!=1: return None
    n=[Fr(str(x)) for x in ns[0]]; return n, sum(a*b for a,b in zip(n,pts[0]))
facets=[]
for S in itertools.combinations(range(9),4):
    h=hyperplane([R[i] for i in S])
    if h is None: continue
    n,b=h; vals=[sum(a*x for a,x in zip(n,R[k]))-b for k in range(9) if k not in S]
    if all(v<0 for v in vals) or all(v>0 for v in vals): facets.append(S)
    assert all(v!=0 for v in vals), 'not general position'
E_exact=set()
for S in facets:
    for i,j in itertools.combinations(S,2): E_exact.add((i,j))
verts=sorted(set(i for S in facets for i in S))
def hi_edge(V,i,j):
    V=np.array([[float(x) for x in v] for v in V]); n,d=V.shape; others=[k for k in range(n) if k not in (i,j)]
    A=np.array([V[k]-V[i] for k in others]); r=linprog(np.zeros(d),A_ub=A,b_ub=-np.ones(len(others)),A_eq=np.array([V[i]-V[j]]),b_eq=np.zeros(1),bounds=[(-1e6,1e6)]*d,method='highs')
    return r.status==0
E_hi=set((i,j) for i,j in itertools.combinations(range(9),2) if hi_edge(R,i,j))
print('facets',len(facets),'vertices',len(verts),'exact edges',len(E_exact),'HiGHS edges',len(E_hi),'agree',E_exact==E_hi)
print('only exact',sorted(E_exact-E_hi),'only HiGHS',sorted(E_hi-E_exact))
