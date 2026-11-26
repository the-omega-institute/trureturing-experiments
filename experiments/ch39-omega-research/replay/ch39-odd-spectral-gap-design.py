"""Support limit for odd sample counts and monotonicity under momentum-grid design."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx

sharp = run_path(str(Path(__file__).with_name("ch39-sharp-hidden-spectral-gap.py")))
ball, roots = sharp["ball"], sharp["gaussian_nodes"]


def radau_spectrum(xs, left):
    order = len(xs)//2
    moment_matrix, denominator = sharp["divided_moments"](xs)
    m = moment_matrix.entries()
    weighted = [m[i+1]-left*m[i] for i in range(2*order)]
    nodes = [left, *roots(weighted,order)]
    assert all(a<b for a,b in zip(nodes,nodes[1:]))
    system = arb_mat([[s**i for s in nodes] for i in range(order+1)])
    weights = system.solve(arb_mat([[v] for v in m[:order+1]])).entries()
    assert min(weights)>0
    for i,actual in enumerate(m):
        residual = sum(w*s**i for w,s in zip(weights,nodes,strict=True))-actual
        assert residual.contains(0) and abs(residual)<arb("1e-20")
    restored = [w*denominator(s) for w,s in zip(weights,nodes,strict=True)]
    for x in xs:
        response=sum(w*ball(x-1)/(1+ball(x-1)*s) for s,w in zip(nodes,restored,strict=True))
        assert (response-sharp["rational_log"](x)).contains(0)
    return (1/nodes[-1]-1).sqrt(), restored[0]


def design_limit(count):
    order = count//2
    # The limiting divided measure is ds, or s ds for odd counts.
    shift = count%2
    moments = [ball(Fraction(1,j+1+shift)) for j in range(2*order)]
    return (1/roots(moments,order)[-1]-1).sqrt()


def verify():
    with ctx.workprec(400):
        for count in (3,5,7):
            xs = [1+Fraction(2*i,count) for i in range(1,count+1)]
            limiting, contact = radau_spectrum(xs,arb(0))
            finite, _ = radau_spectrum(xs,ball(Fraction(1,1024)))
            tighter, _ = radau_spectrum(xs,ball(Fraction(1,2048)))
            assert 0<finite<tighter<limiting and contact>0
            print(count,"critical M/M_N",limiting.str(15),"nonzero contact weight",contact.str(12),
                  "finite-atom M/M_N at 1/1024",finite.str(14),"at 1/2048",tighter.str(14))
        for count in (2,3,4,5,6,7,8):
            lower = design_limit(count)
            xs = [1+Fraction(i,16*count) for i in range(1,count+1)]
            if count%2:
                near, _ = radau_spectrum(xs,arb(0))
            else:
                mm,_ = sharp["divided_moments"](xs)
                near = (1/roots(mm.entries(),count//2)[-1]-1).sqrt()
            assert near>lower
            print(count,"momentum design infimum",lower.str(14),"narrow grid",near.str(14))


if __name__ == "__main__":
    verify()
