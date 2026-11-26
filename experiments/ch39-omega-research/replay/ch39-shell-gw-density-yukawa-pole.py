"""Uniformly attainable equal-site light data and a matrix-flow pole on two RH phase windows."""
from pathlib import Path
from runpy import run_path
import numpy as np
from flint import acb_mat, arb

folder = Path(__file__).resolve().parent
flow = run_path(str(folder/"ch39-shell-gw-first-heavy-flow.py"))
density = run_path(str(folder/"ch39-shell-gw-positive-density-cp.py"))
fit = flow["fit"]
rows,count,hs,adj = (fit[k] for k in ("rows","count","hs","adj"))


def equal_fit(P,up,down,e):
    old = fit["common"]["old"]
    U,A,B = [fit["orbit"]["exact_frame"](np.array(z).reshape(3,3)) for z in old["units"]]
    target = acb_mat([[z] for z in (U*up*A).entries()+(U*down*B).conjugate().entries()])
    target *= arb("12.5")/arb("7.9")
    assert count(rows(P*adj(P)),arb(1)/400)==0
    T = adj(P)*(P*adj(P)).inv()
    z = density["cp"]["trial"](P,target)
    x = z+T*(target-P*z)
    N = old["m"]["eye"](54)-T*P
    derivative = acb_mat([[2*x[j,0].conjugate() if j//9==i else 0 for j in range(54)]
                          for i in range(6)])
    G = derivative*N*adj(derivative)
    real = acb_mat([[z.real for z in row] for row in rows(G)])
    C = arb(".45")
    assert count(rows(real),1/C**2)==0
    masses = [sum(abs(x[9*i+j,0])**2 for j in range(9)) for i in range(6)]
    residual = sum(abs(t-arb("12.5")**2/6).abs_upper()**2 for t in masses).sqrt()
    shift = 20*e/(1-20*e)
    error = shift*hs(x,x).real.sqrt()
    eta = residual+2*max(masses).sqrt()*error+error**2
    slope = 2*max(masses).sqrt()*shift*C+2*error*(1+shift)*C
    rate = slope+4*((1+shift)*C)**2*eta
    assert rate < arb(".5")
    return rate


def check(index,width):
    e = density["window"](index,width)
    Q,m = flow["masses"](index)
    # The earlier certified Hermitian perturbation of each gap is <= 4*width.
    # The same total gauge exposure pins the common M between the two
    # reciprocal square-root scalings; the common light offset is positive.
    gap_lower = (m[0]**2-125**2)/arb("1.2e12")**2
    eps = 4*width/gap_lower
    factor = ((1+eps)/(1-eps)).sqrt()
    stretch = (1/factor).union(factor)
    first,last = m[0]*stretch,m[-1]*stretch
    assert last < Q and first>1000
    end = (first/1000).log()/flow["loop"]
    up,down,_ = flow["evolve"](end)
    F,D = fit["readouts"](index)
    P = acb_mat(rows(F)+rows(D))
    rate = equal_fit(P,up,down,e)
    b = arb("5.14" if index==43 else "5.06")
    u,d = (hs(y,y).real.sqrt() for y in (up,down))
    limit = arb(".01080")
    assert count(rows(adj(F)*F-60*adj(D)*D),limit)==54
    assert u>e*b and (u-e*b)**2-60*(d+e*b)**2>limit*b**2
    norms = [hs(Y,Y).real.sqrt()-3*e for Y in (F,D)]
    assert min(norms)>0
    gfirst = flow["gauge"](end,(flow["start"],flow["beta"]))
    unified = flow["gauge"]((Q/1000).log()/flow["loop"],(flow["start"],flow["beta"]))[2]
    caps = [arb.max(v,unified) for v in gfirst]
    assert all(2*sum(z*v for z,v in zip(row,caps))<9 for row in flow["C"][:2])
    S = b**2*sum(t**2 for t in norms)*(-9*(last/first).log()/flow["loop"]).exp()
    assert S > 9
    pole = last*(flow["loop"]/9*(S/(S-9)).log()).exp()
    assert pole < Q
    print("zero",index,"square half-width",width.str(2),"equal-site light fit contraction <",rate.upper().str(6),
          "unfiltered pole by",pole.upper().str(9),"GeV")


if __name__=="__main__":
    for n,w in ((43,"3e-8"),(195,"1e-8")):
        check(n,arb(w))
