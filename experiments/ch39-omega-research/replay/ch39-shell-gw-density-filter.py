"""Simultaneous local filtered matching on two positive-density RH zero windows."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb_mat, arb

w = Path(__file__).resolve().parent
inert = run_path(str(w/"ch39-shell-gw-two-inert-heavy.py"))
local = run_path(str(w/"ch39-shell-gw-local-yukawa-filter.py"))
window = run_path(str(w/"ch39-shell-gw-density-yukawa-pole.py"))
flow, fit = inert["flow"],inert["fit"]
rows, count, adj, hs = (fit[k] for k in ("rows","count","adj","hs"))


def frob_upper(M):
    return sum(abs(z).abs_upper()**2 for z in M.entries()).sqrt()


def triple(index,width,e,up,down,F,D):
    P = acb_mat(rows(F)+rows(D))
    assert count(rows(P*adj(P)),arb(1)/400)==0
    T = adj(P)*(P*adj(P)).inv()
    N = fit["common"]["old"]["m"]["eye"](54)-T*P
    old = fit["common"]["old"]
    U,A,B = [fit["orbit"]["exact_frame"](np.array(z).reshape(3,3)) for z in old["units"]]
    target = acb_mat([[z/inert["g"],0,0] for z in (U*up*A).entries()+(U*down*B).conjugate().entries()])
    grams = inert["fingerprints"](index)
    trial = inert["trial"](P,target,N,grams)
    X = trial+T*(target-P*trial)
    initial2 = arb(0)
    gradients = []
    for j in range(6):
        site = acb_mat(rows(X)[9*j:9*j+9])
        gram = adj(site)*site
        assert count(rows(gram),arb(".42")**2)==3
        initial2 += frob_upper(gram-grams[j])**2
        for axis in inert["axes"]():
            embed = (2*site*axis).table()
            block = acb_mat([embed[i%9] if i//9==j else [0,0,0] for i in range(54)])
            gradients.append((N*block).entries())
    G = acb_mat(gradients)
    real = acb_mat([[z.real for z in row] for row in rows(G*adj(G))])
    sigma = arb(".12" if index==43 else ".15")
    assert count(rows(real),sigma**2)==0
    # Choose an orthonormal frame in the whole perturbed two-mode band using
    # the polar factor of the projected central frame.
    delta = 4*width
    basis_error = 3*delta
    shift = 20*e/(1-20*e)
    move = shift*frob_upper(X)
    eta = initial2.sqrt()+2*basis_error+basis_error**2+arb(".84")*move+move**2
    sigma1 = sigma-arb(".84")*(2*shift)-2*move
    assert sigma1>0
    rate = 4*eta/sigma1**2
    # The map on a ball of radius 2*eta has image norm <= eta + rate*eta.
    assert rate < 1
    return rate


def check(index,width):
    Q,m = flow["masses"](index)
    model = flow["full"]["model"]
    if index==43:
        s,f,_,_=model["m"]["transported_spectra"]()
    else:
        s,f,_=model["scalar_models"]()
    assert count(rows(f),arb(".04"))==count(rows(f),arb("1.05"))==2
    gaps = flow["full"]["all_levels"](s,f)
    delta = 4*width
    assert gaps[1]+delta < 1 < gaps[2]-delta and gaps[-1]+delta < 7
    # Individual central eigenvectors have residual <=delta in the perturbed
    # operator; the separation from the exterior of the cluster exceeds one.
    assert arb("1.05")-arb(".04")-delta > 1
    eps = delta/gaps[0]
    stretch = ((1+eps)/(1-eps)).sqrt()
    band = (1/stretch).union(stretch)
    first,third = m[0]*band,m[2]*band
    end = (first/1000).log()/flow["loop"]
    data = flow["evolve"](end)
    F,D = fit["readouts"](index)
    e = window["density"]["window"](index,width)
    rate = triple(index,width,e,*data[:2],F,D)
    # Every mode above the protected pair lies in the continuous stop band.
    gain = 1/local["T4"](arb(4)/3)
    total = sum((hs(A,A).real.sqrt()+3*e)**2 for A in (F,D))
    y = sum(hs(M,M).real for M in data)
    early = (third/first).log()/flow["loop"]
    assert 10*early*y<1
    y = y/(1-10*early*y)+inert["g"]**2*total*gain**2
    late = (Q/third).log()/flow["loop"]
    assert 10*late*y<1
    y /= 1-10*late*y
    assert y < arb("2.59")
    print("zero",index,"whole-window triple contraction below",rate.upper().str(8),
          "whole-window running Yukawa squared norm below",y.upper().str(8))


if __name__=="__main__":
    for n,width in ((43,"3e-8"),(195,"1e-8")):
        check(n,arb(width))
