"""Neutral states forced by the adjoint lift and commuting vertex Higgs backgrounds."""
from pathlib import Path
import runpy
from flint import acb,acb_mat,arb
import numpy as np

lift = runpy.run_path(str(Path(__file__).with_name("ch39-su3-root-clock-lift.py")))
short,m,shell,q = (lift[k] for k in ("short","m","sh","q"))
I = lift["I"]
sigma = [acb_mat([[0,1],[1,0]]), acb_mat([[0,acb(0,-1)],[acb(0,1),0]]),
         acb_mat([[1,0],[0,-1]])]
def adjoint(g):
    gd = g.conjugate().transpose()
    result = []
    for a in sigma:
        row = []
        for b in sigma:
            v = a*g*b*gd
            z = (v[0,0]+v[1,1])/2
            assert z.imag.contains(0)
            row.append(z.real)
        result.append(row)
    return result

turn = {d:adjoint(U) for d,U in ((724,lift["U"]),(804,lift["V"]),(0,I))}
weights = short["base"]["weights"]
def matrix(weights,dim):
    H = acb_mat(dim*19,dim*19)
    ids = m["source_ids"]
    for i,j,d,_ in shell["edges"]:
        w = 1+weights[ids[i,j]] if (i,j) in ids else q(1)
        for a in range(dim):
            for v in (i,j):
                H[dim*v+a,dim*v+a] += w
            for b in range(dim):
                R = q(1) if dim==1 else turn.get(d,turn[0])[a][b]
                H[dim*i+a,dim*j+b] = -w*R
                H[dim*j+b,dim*i+a] = -w*R
    return H

def cutoff(H,j,lower):
    n = H.nrows()
    mid = np.array([[float(H[a,b].real.mid()) for b in range(n)] for a in range(n)])
    ev = np.linalg.eigvalsh(mid)[j]
    bound = arb(str(ev))+q("-.001" if lower else ".001")
    count = lift["source"]["inertia"](H,bound)
    if lower:
        assert count<=j
    else:
        assert count>j
    return bound

def check():
    d0 = short["base"]["data"]["d"]
    eta = q("1.018")**2
    assert d0/eta>short["r"]**2
    _,cap = short["base"]["arc_budget"]()
    k = arb.pi()*(cap*eta/d0).sqrt()
    slope = 18*k*(k/2).tan()/cap
    for dim,j in ((3,0),(1,1)):
        L = matrix([v.lower() for v in weights],dim)
        R = matrix([v.upper() for v in weights],dim)
        lo = cutoff(L,j,True)/(1+slope)
        hi = cutoff(R,j,False)
        top = cutoff(R,19*dim-1,False)
        assert top<cap
        print("neutral block",dim,"first positive mass squared",lo.lower().str(12),hi.upper().str(12),
              "all vertex levels below",top.upper().str(12))
    assert cap<200
    assert min((v+p).lower() for v,p in zip(short["potential"],short["shiftbox"]))>200
    print("exact neutral levels below the charged floor at 200 M^2:",19*(3+1))
if __name__=="__main__":
    check()
