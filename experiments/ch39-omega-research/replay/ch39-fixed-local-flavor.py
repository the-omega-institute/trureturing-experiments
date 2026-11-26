"""Reuse one dyadic four-parameter action throughout the gauge candidate window."""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb

arc = run_path(str(Path(__file__).with_name("ch39-local-flavor-arc-count.py")))
native = arc["native"]
shell = arc["shell"]
basis = arc["basis"]
identity = acb_mat([[1,0],[0,1]])
radial = acb_mat([[int(i == j) for j in range(3)] for i in range(3)]) - basis*basis.transpose()
arcs = arc["arc_projectors"]()
codes = (1486,5202,5069,6322)


def mixing(t, coefficients):
    g,h,k,phase = coefficients
    action = g*t+h*arcs[32]+k*arcs[152]+phase*identity
    mean = action.trace().real/2
    reduced = action-mean*identity
    gap = (reduced[0,0].real**2+abs(reduced[0,1])**2).sqrt()
    assert gap > 0
    u = shell["cis"](mean)*(gap.cos()*identity + acb(0,gap.sin()/gap)*reduced)
    return radial+basis*u*basis.transpose()


def errors(w):
    u = abs(w[0,2])
    c = (1-u*u).sqrt()
    s,t = abs(w[0,1])/c, abs(w[1,2])/c
    a,b = (1-s*s).sqrt(),(1-t*t).sqrt()
    j = (w[0,0]*w[1,0].conjugate()*w[1,1]*w[0,1].conjugate()).imag
    sd = j/(a*b*c**2*s*t*u)
    cd = (abs(w[1,0])**2-(s*b)**2-(a*t*u)**2)/(2*a*b*s*t*u)
    angle = acb(cd,sd).arg()
    return [abs(x-y) for x,y in zip((s,t,u,angle),shell["centers"])]


def verify():
    coeff = [arb(n)/2048 for n in codes]
    varying = [v+arb(0,arb(2)**-28) for v in coeff]
    tolerances = [shell["sine_radius"]]*3+[shell["phase_radius"]]
    winners = []
    runner = []
    for idx,z in enumerate(acb.zeta_zeros(10815,1278),10815):
        response = native["plane_response"](z.imag)
        vals = errors(mixing(response,varying))
        if all(e < t for e,t in zip(vals,tolerances)):
            winners.append(idx)
            center = errors(mixing(response,coeff))
            print(idx,"CKM deviations",*(v.str(7) for v in center))
        else:
            assert any(e > t for e,t in zip(vals,tolerances))
            runner.append((max((e/t).lower() for e,t in zip(vals,tolerances)),idx))
    print("fixed action, matching indices",winners)
    assert winners == [11431]
    assert all(lower > arb("10.1272") for lower,_ in runner)
    print("smallest certified competing scaled-error lower bound",min(runner))


if __name__ == "__main__":
    verify()
