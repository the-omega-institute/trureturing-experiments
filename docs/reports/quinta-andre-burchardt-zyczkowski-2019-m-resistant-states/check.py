# Check of the flag-coded construction for Quinta–André–Burchardt–Życzkowski's m-resistance conjecture
# (Zhang–Han–Shi–Zhang arXiv:2505.06567 Conjecture 1): for N >= 3 and 0 <= m <= N-2, with k = N - m, s = C(N,k), d = 2s,
# psi = (2s)^(-1/2) sum_{S in C([N],k)} sum_{b in {0,1}} tensor_i |S, b*[i in S]>.
import itertools, sys
import numpy as np
from fractions import Fraction
def states(N,m):
    k=N-m; E=list(itertools.combinations(range(N),k)); s=len(E)
    idx={(S,b):j for j,(S,b) in enumerate((S,b) for S in E for b in (0,1))}   # local basis index
    terms=[tuple(idx[(S, b if i in S else 0)] for i in range(N)) for S in E for b in (0,1)]
    return E,s,terms
def reduced(terms,s,keep):
    # rho_K entries: sum over pairs of terms agreeing on lost sites
    lost=[i for i in range(len(terms[0])) if i not in keep]
    rho={}
    amp=Fraction(1,2*s)
    for t in terms:
        for u in terms:
            if all(t[i]==u[i] for i in lost):
                key=(tuple(t[i] for i in keep),tuple(u[i] for i in keep))
                rho[key]=rho.get(key,0)+amp
    return rho
def is_diagonal(rho): return all(a==b for (a,b) in rho)
def pt_min_eig(rho,keep_len,site):
    # partial transpose at position `site` of the kept tuple; restrict to the support basis to get a finite matrix
    pt={}
    for (a,b),v in rho.items():
        a2=list(a); b2=list(b); a2[site],b2[site]=b[site],a[site]
        pt[(tuple(a2),tuple(b2))]=pt.get((tuple(a2),tuple(b2)),0)+v
    basis=sorted({x for k in pt for x in k}); ix={x:i for i,x in enumerate(basis)}
    M=np.zeros((len(basis),len(basis)))
    for (a,b),v in pt.items(): M[ix[a],ix[b]]+=float(v)
    return np.linalg.eigvalsh(M).min()
ok=True
for N in range(3,7):
    for m in range(0,N-1):
        E,s,terms=states(N,m)
        assert len(set(terms))==2*s   # orthonormal summands
        worst_pt=0; diag_ok=True
        for L in itertools.combinations(range(N),m):
            keep=[i for i in range(N) if i not in L]
            rho=reduced(terms,s,keep)
            worst_pt=min(worst_pt,min(pt_min_eig(rho,len(keep),j) for j in range(len(keep))))
        for L in itertools.combinations(range(N),m+1):
            keep=[i for i in range(N) if i not in L]
            rho=reduced(terms,s,keep); diag_ok&=is_diagonal(rho)
        good = worst_pt < -1e-12 and diag_ok
        ok&=good
        print(f'N={N} m={m} d={2*s}: every |L|=m reduced state has a negative partial transpose (min eig {worst_pt:.4g}); every |L|=m+1 reduced state is diagonal in a product basis: {diag_ok}', 'OK' if good else 'FAIL')
print('ALL_OK' if ok else 'FAIL'); sys.exit(0 if ok else 1)
