#!/usr/bin/env python3
"""Zhang-Xie (arXiv:2601.01858v2), open question after Corollary 5.5: classify the n-tuples of unit vectors in C^2 with
Im tr(psi_1 ... psi_n) = tau_n = cos^n(pi/n) sec^(n-1)(pi/(2(n-1))). Proposed answer: the unitary-and-phase orbit of
v_j = (sqrt(1-t_n), sqrt(t_n) e^(2 pi i j/n)), t_n = (1 - tan(pi/(2(n-1)))/tan(pi/n))/2.
For n = 3..8: (1) the family attains tau_n; (2) 40 random-restart local maximisations of Im tr reach tau_n and never
exceed it beyond 1e-12; (3) the best tuple found matches the family's |Gram| matrix up to cyclic relabelling (1e-6).
Exit 0 and RESULT PASS when all hold.
"""
import sys
import numpy as np
rng=np.random.default_rng(1858)
OK=[]
def imdelta(V):  # V: n x 2 unit rows; tr(psi_1...psi_n) = prod_j <psi_j, psi_{j+1}>
    n=len(V); p=1+0j
    for j in range(n): p*=np.vdot(V[j],V[(j+1)%n])
    return p.imag
def tau(n): return np.cos(np.pi/n)**n/np.cos(np.pi/(2*(n-1)))**(n-1)
def family(n):
    a=np.pi/n; b=np.pi/(2*(n-1)); t=(1-np.tan(b)/np.tan(a))/2; w=np.exp(2j*np.pi/n)
    return np.array([[np.sqrt(1-t),np.sqrt(t)*w**j] for j in range(n)])
def unit(x):
    V=(x[:,:2]+1j*x[:,2:]); return V/np.linalg.norm(V,axis=1,keepdims=True)
for n in range(3,9):
    fam=imdelta(family(n)); best=-1; bestV=None
    for trial in range(40):
        x=rng.normal(size=(n,4)); g=imdelta(unit(x)); step=0.5
        for it in range(4000):
            y=x+rng.normal(size=(n,4))*step; gy=imdelta(unit(y))
            if gy>g: x,g=y,gy; step*=1.2
            else: step*=0.97
            if step<1e-10: break
        if g>best: best,bestV=g,unit(x)
    # orbit invariant: all |<psi_j,psi_{j+1}>| equal and the Gram matrix matches the family's up to phases
    def gram_abs(V): return np.abs(V.conj()@V.T)
    G1=gram_abs(bestV); G0=gram_abs(family(n))
    # compare up to cyclic shift / reflection of the index
    dists=[]
    for s in range(n):
        for rev in (False,True):
            idx=[(s-j)%n if rev else (s+j)%n for j in range(n)]
            dists.append(np.abs(G1[np.ix_(idx,idx)]-G0).max())
    ok_n = abs(fam-tau(n))<1e-12 and abs(best-tau(n))<1e-9 and best-tau(n)<1e-12 and min(dists)<1e-6
    OK.append(ok_n)
    print(f"ok={ok_n}", end=" "); print(f'n={n} tau={tau(n):.12f} family={fam:.12f} max_found={best:.12f} (max-tau={best-tau(n):.2e}) |Gram| distance to family orbit={min(dists):.2e}')
print("RESULT PASS" if all(OK) else "RESULT FAIL"); sys.exit(0 if all(OK) else 1)
