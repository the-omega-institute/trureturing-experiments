# Brute force for Lillo-Rosas arXiv:2603.28905: theta_n (principal filters [P,1] iso to some Pi(m), m<=n)
# and gamma_n (principal ideals [0,P] iso to some Pi(m), 1<=m<=n). Priority forests on labels 0..n:
# consecutive increasing component trees; order by edge inclusion, plus a top element.
import itertools, sys, networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher
from math import factorial
def compositions(N):
    for k in range(1,N+1):
        for cuts in itertools.combinations(range(1,N),k-1):
            b=(0,)+cuts+(N,); yield [b[i+1]-b[i] for i in range(k)]
def forests(n):
    N=n+1; out=[]
    for comp in compositions(N):
        starts=[sum(comp[:i]) for i in range(len(comp))]
        choices=[]
        for s,L in zip(starts,comp):
            per=[[ (p,v) for p in range(s,v)] for v in range(s+1,s+L)]
            choices.append(list(itertools.product(*per)) if per else [()])
        for pick in itertools.product(*choices):
            out.append(frozenset(e for t in pick for e in t))
    return out
def poset(n):
    F=forests(n); els=F+['TOP']
    return els
def hasse(elements):
    # elements: list of frozensets or 'TOP'; order by inclusion, TOP above all
    G=nx.DiGraph(); idx={e:i for i,e in enumerate(elements)}
    G.add_nodes_from(range(len(elements)))
    fs=[e for e in elements if e!='TOP']
    for a in fs:
        for b in fs:
            if len(b)==len(a)+1 and a<b: G.add_edge(idx[a],idx[b])
    if 'TOP' in idx:

        for a in fs:
            if not any(a<b for b in fs): G.add_edge(idx[a],idx['TOP'])
    return G
cache={}
def Pi(m):
    if m not in cache: cache[m]=hasse(poset(m))
    return cache[m]
def iso(G,H):
    if G.number_of_nodes()!=H.number_of_nodes() or G.number_of_edges()!=H.number_of_edges(): return False
    return DiGraphMatcher(G,H).is_isomorphic()
NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 5
for n in range(1,NMAX+1):
    F=forests(n)
    theta=1  # P = TOP: [TOP,TOP] is one element, not iso to any Pi(m) (Pi(0) has 2 elements) -> excluded below
    theta=0; gamma=0; gamma0=0
    for P in F+['TOP']:
        if P=='TOP':
            up=['TOP']; down=F+['TOP']
        else:
            up=[Q for Q in F if P<=Q]+['TOP']; down=[Q for Q in F if Q<=P]
        Gu=hasse(up); Gd=hasse(down) if P!='TOP' else Pi(n)
        if any(iso(Gu,Pi(m)) for m in range(0,n+1)): theta+=1
        if any(iso(Gd,Pi(m)) for m in range(1,n+1)): gamma+=1
        if any(iso(Gd,Pi(m)) for m in range(0,n+1)): gamma0+=1
    left=sum(factorial(k) for k in range(n+1))
    print(n,len(F)+1,"theta",theta,"left_factorial(n+1)",left,"gamma(m>=1)",gamma,"gamma(m>=0)",gamma0,"n^2-n+2",n*n-n+2,flush=True)
