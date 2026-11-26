"""Certify a common fitted doublet profile with equal norm at all six shell sites."""
from pathlib import Path
from runpy import run_path
import numpy as np
from flint import acb, acb_mat, arb

common=run_path(str(Path(__file__).with_name("ch39-shell-gw-195-common-field.py")))
model=common["flavor"]
adj,hs,rows=common["adj"],common["hs"],common["rows"]


def numeric_trial(P,target):
    A=model["array"](P)
    T=A.conj().T@np.linalg.inv(A@A.conj().T)
    x=T@model["array"](target).ravel()
    kernel=np.eye(54)-T@A
    for _ in range(18):
        sites=x.reshape(6,9)
        sq=np.sum(abs(sites)**2,axis=1)
        residual=sq[:-1]-sq[-1]
        if np.max(abs(residual))<1e-12:
            break
        D=np.zeros((5,54),complex)
        for j in range(5):
            D[j,9*j:9*(j+1)]=2*sites[j].conj()
            D[j,-9:]=-2*sites[-1].conj()
        D=D@kernel
        x-=kernel@D.conj().T@np.linalg.solve((D@D.conj().T).real,residual)
    return acb_mat([[acb(float(v.real),float(v.imag))] for v in x])


def certify(n):
    F=common["old"]["F"] if n==43 else common["operators"]()[1]
    t=common["native"]["relative_phase"](arb(676))
    P=acb_mat(rows(F)+rows(common["conjugate"](F,t)))
    target=common["arc"]["central_target"]()
    seed=numeric_trial(P,target)
    T=adj(P)*(P*adj(P)).inv()
    x=seed+T*(target-P*seed)
    kernel=model["eye"](54)-T*P
    blocks=[acb_mat([[x[9*i+j,0]] for j in range(9)]) for i in range(6)]
    sq=[hs(z,z).real for z in blocks]
    residual=sum((a-sq[-1])**2 for a in sq[:-1]).sqrt()
    D=acb_mat([[2*x[j,0].conjugate()*(int(j//9==i)-int(j//9==5))
                for j in range(54)] for i in range(5)])
    V=kernel*adj(D)
    gram=D*V
    real=acb_mat([[gram[i,j].real for j in range(5)] for i in range(5)])
    C=V*real.inv()
    norm=hs(C,C).real.sqrt()
    radius=2*residual
    assert 2*arb(5).sqrt()*norm**2*radius<arb(".5")
    upper=hs(x,x).real.sqrt()+norm*radius
    assert upper<arb("12.29" if n==43 else "12.34")<4*arb.pi()
    print("zero",n,"six equally populated sites, exact joint fit with strength below",upper.upper().str(10))


if __name__=="__main__":
    for n in (43,195):
        certify(n)
