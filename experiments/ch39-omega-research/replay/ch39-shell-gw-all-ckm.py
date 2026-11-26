"""Verify all down-sector left Gram matrices in a ball, at fixed up matrix and site norms."""
from pathlib import Path
from runpy import run_path
import numpy as np
from flint import acb, acb_mat, arb

density=run_path(str(Path(__file__).with_name("ch39-shell-gw-positive-density-cp.py")))
cp=density["cp"]
common,model,adj,hs,rows=(cp[key] for key in ("common","model","adj","hs","rows"))
count=common["count"]


def real_stack(A):
    return acb_mat([[z.real for z in row] for row in rows(A)]+
                   [[z.imag for z in row] for row in rows(A)])


def real_action(A):
    return acb_mat([[z.real for z in row]+[-z.imag for z in row] for row in rows(A)]+
                   [[z.imag for z in row]+[z.real for z in row] for row in rows(A)])


def norm_below(A,bound):
    gram=adj(A)*A
    assert count(rows(gram),bound**2)==gram.nrows()


def subspace(lift):
    # Binary proposals for seven complex output directions. The Gram and all
    # lift norms are certified after converting the coefficients to exact reals.
    vectors=np.linalg.svd(model["array"](lift).real,full_matrices=False)[2][:2]
    constraints=vectors[:,:9]-1j*vectors[:,9:]
    basis=np.linalg.svd(constraints,full_matrices=True)[2][2:].conj().T
    H=acb_mat([[acb(float(z.real),float(z.imag)) for z in row] for row in basis])
    assert count(rows(adj(H)*H),arb(".9999")**2)==0
    return H


def certify(index,width):
    op=density["window"](index,width)
    F=common["old"]["F"] if index==43 else common["operators"]()[1]
    t=common["native"]["relative_phase"](arb(676))
    P=acb_mat(rows(F)+rows(common["conjugate"](F,t)))
    y=common["arc"]["central_target"]()
    target=acb_mat([[y[i,0] if i<9 else 0] for i in range(18)])
    seed=cp["trial"](P,target)
    T=adj(P)*(P*adj(P)).inv()
    assert count(rows(P*adj(P)),arb(1)/400)==0
    x=seed+T*(target-P*seed)
    N=model["eye"](54)-T*P
    D=acb_mat([[2*x[j,0].conjugate() if j//9==i else 0 for j in range(54)] for i in range(6)])
    Z=N*adj(D)
    gram=D*Z
    real=acb_mat([[gram[i,j].real for j in range(6)] for i in range(6)])
    C=Z*real.inv()
    c=arb(".60")
    norm_below(real_stack(C),c)
    Td=acb_mat([row[9:] for row in rows(T)])
    lift=real_action(Td)-real_stack(C)*acb_mat(rows(real_action(D*Td))[:6])
    H=subspace(lift)
    L=lift*real_action(H)
    norm_below(L,arb("14.2"))
    for i in range(6):
        idx=list(range(9*i,9*i+9))+list(range(54+9*i,63+9*i))
        block=acb_mat([rows(L)[j] for j in idx])
        norm_below(block,arb("8.1"))
    eps=arb(".01402")/arb(".9999")
    total,site=arb("14.2")*eps,arb("8.1")*eps
    masses=[sum(abs(x[9*i+j,0])**2 for j in range(9)) for i in range(6)]
    amplitude=max(masses).sqrt()
    residual=sum((v-cp["strength"]**2/6)**2 for v in masses).sqrt()
    q=20*op/(1-20*op)
    displacement=q*(hs(x,x).real.sqrt()+total)
    eta=residual+total**2+2*(amplitude+site)*displacement+displacement**2
    slope=2*amplitude*c*q+2*(site+displacement)*c*(1+q)
    radius=2*eta
    contraction=slope+2*(c*(1+q))**2*radius
    assert contraction<arb(".5")
    print("zero",index,"uniform contraction for the entire down Gram ball below",
          contraction.upper().str(8))


if __name__=="__main__":
    for n,w in ((43,"3e-8"),(195,"1e-8")):
        certify(n,arb(w))
