"""Certify the new winner after widening the gauge input box, and test the two port readouts."""
from pathlib import Path
from runpy import run_path
from flint import acb, arb, arb_series, ctx

root = Path(__file__).parent
gauge = run_path(str(root/"ch39-gauge-width-transition.py"))
finite = run_path(str(root/"ch39-finite-arithmetic-height.py"))
q,weight,hardy = gauge["q"],finite["weight"],finite["hardy"]


def verify():
    ctx.prec,ctx.cap = 144,3
    gauge["robust_exclusion"](q(".001030"),q("26.8642"),10567)
    winner = acb.zeta_zero(10854).imag
    best = weight(winner)
    envelope = finite["source"]["modulus"](-arb(59)/108)
    assert envelope/(arb(10613)**2+arb(1)/4)<best
    start = int(arb(10567).zeta_nzeros().unique_fmpz())+1
    last = int(arb(10613).zeta_nzeros().unique_fmpz())
    assert start==10812 and 10854<=last
    for n,z in enumerate(acb.zeta_zeros(start,last-start+1),start):
        assert n==10854 or weight(z.imag)<best

    t = arb_series([winner,1])
    critical = abs(weight(t)[1]/best/t.riemann_siegel_z()[1])
    assert arb(".00548833")<critical<arb(".00548835")
    assert arb(1)/208 < critical < arb(1)/13
    near,eps,slack=arb(1)/1024,arb("3.65e-22"),arb("1e-10")
    for j in range(-32,32):
        cell = winner+near*(2*j+1)/64+arb(0,near/64)
        v = arb_series([cell,1])
        assert abs(hardy(cell))+eps < arb(".02")
        assert v.riemann_siegel_z()[1] < -7
        assert abs(weight(v)[1])/best < arb(".06")
    stack = [(arb(10567),(winner-near).upper()),((winner+near).lower(),arb(10613))]
    while stack:
        lo,hi = stack.pop()
        cell = lo.union(hi)
        gain = weight(cell)/best
        if gain<1-slack:
            continue
        z = abs(hardy(cell))+arb(0,eps)
        if gain<(1-slack)*(1+z/(13+12*z)):
            continue
        assert hi-lo>arb("1e-9")
        middle = (lo+hi)/2
        stack.extend([(lo,middle),(middle,hi)])
    assert envelope/(arb(10613)**2+arb(1)/4)<best*(1-slack)

    grid = 2**64
    center = int((winner*grid).floor().unique_fmpz())
    scores = {}
    for m in range(center-6,center+8):
        g = arb(m)/grid
        z = abs(finite["finite_value"](g))
        scores[m] = weight(g)/(1+z/(13+12*z))
    top = max(scores,key=lambda m:scores[m].mid())
    assert all(scores[top]>s for m,s in scores.items() if m!=top)
    level = scores[top].lower()/best
    coeff = 1/(13+12*arb(".02"))
    radius = (1-level+level*coeff*eps)/(7*level*coeff-arb(".06"))
    assert 0<radius<6/arb(grid)
    print("1011 through 1030 ppm: zero 10854, local critical strength",critical.str(14),"binary winner",top,
          "offset in grid units", (top-winner*grid).str(8))


if __name__ == "__main__":
    verify()
