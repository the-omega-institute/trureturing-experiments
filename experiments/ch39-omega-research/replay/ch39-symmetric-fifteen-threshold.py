"""对称十五维物质塔的正阈值、最小正向表示和所需份数。"""
from collections import defaultdict
from fractions import Fraction
from itertools import product
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_poly, ctx

matter = run_path(str(Path(__file__).with_name("ch39-charge-matter-determinant-order.py")))
soft = matter["soft"]


def small_representations():
    # Weyl dimension for SU(5) is increasing in each Dynkin label.
    def dimension(labels):
        return prod(Fraction(j-i+sum(labels[i:j]), j-i)
                    for i in range(4) for j in range(i+1, 5))
    assert all(dimension(tuple(3*int(i == j) for i in range(4))) > 15 for j in range(4))
    found = sorted(dimension(x) for x in product(range(3), repeat=4) if dimension(x) <= 15)
    assert found == [1, 5, 5, 10, 10, 15, 15]


def derivative_bounds(tables):
    total = defaultdict(arb)
    for a, b, c, factor in ((4, 6, 1, 7), (1, 4, 6, -4)):
        delta = {p: value-tables[b][p] for p, value in tables[a].items() if matter["support"](*p)}
        w = matter["mass_wronskian"](delta, tables[b])
        for (i, j), value in w.items():
            for (l, k), coeff in tables[c].items():
                if (l, k) != (0, 0):
                    total[i+l, j+k] += factor*value*coeff
    assert len(total) == 1585 and all(value > 0 for value in total.values())
    # At x=0 all three determinants have a factor e which cancels in the ratio.
    a, b, c = [arb_poly([tables[q][i, 0] for i in range(1, 19)]) for q in (1, 4, 6)]
    derivative = 4*a.derivative()*b*c+7*c.derivative()*a*b-11*b.derivative()*a*c
    # The degree-50 term is zero since 4+7=11 and all degrees are 17.
    assert derivative[50].contains(0) and all(derivative[k] > 0 for k in range(50))
    return (11*b[0].log()-4*a[0].log()-7*c[0].log())/(5*arb.pi())


def verify():
    ctx.prec = 224
    small_representations()
    height = acb.zeta_zero(11431).imag
    base, neutral, charged = soft["vector_pencil"](height)
    blocks = {q:matter["nonspokes"](q, height) for q in (1, 4, 6)}
    tables = {q:matter["coefficients"](base, b) for q, b in blocks.items()}
    ceiling = derivative_bounds(tables)
    scale = (arb("1e17")/arb("91.1876")).log()
    em, weak = arb("127.930"), arb(".23129")
    sm = [3*(1-weak)*em/5-arb("4.1")*scale/(2*arb.pi()),
          weak*em+arb(19)*scale/(12*arb.pi()), 1/arb(".1180")+7*scale/(2*arb.pi())]
    contrast = sm[0]-3*sm[1]+2*sm[2]
    radius = matter["geo"]["theta"](804)/matter["geo"]["theta"](12)
    for weight, count in ((arb(1)/10000, 37), (arb(1), 63)):
        def logs(x):
            return {q:(base+weight*b+x*soft["identity"]).det().real.log() for q, b in blocks.items()}
        def shift(x):
            d = logs(x)
            return (11*d[4]-4*d[1]-7*d[6])/(5*arb.pi())
        bound = ceiling if count == 37 else shift(arb(0))
        assert contrast+(count-1)*bound < 0
        lo, hi = arb(0), arb(1)
        assert contrast+count*shift(lo) > 0 and contrast+count*shift(hi) < 0
        for _ in range(80):
            mid = (lo+hi)/2
            value = contrast+count*shift(mid)
            if value > 0:
                lo = mid
            else:
                assert value < 0
                hi = mid
        mass2 = lo.union(hi)
        exposures = {q:19*radius.log()-d/2 for q, d in logs(mass2).items()}
        lift = acb_mat([[arb(1)/19]*19 for _ in range(19)])
        v0 = 18*radius.log()-(base+weight*neutral+lift).det().real.log()/2
        v5 = 19*radius.log()-(base+weight*charged).det().real.log()/2
        indices = ((arb(1)/10, arb(8)/5, arb(9)/5), (arb(3)/2, 0, 2), (1, arb(5)/2, 0))
        absolute = [sm[i]+7*((0,2,3)[i]*v0+(5,3,2)[i]*v5)/(4*arb.pi())
                    -2*count*sum(t*exposures[q] for t, q in zip(row, (1,4,6)))/(3*arb.pi())
                    for i, row in enumerate(indices)]
        assert all(value < -678 for value in absolute)
        print("Count, fitted bare mass squared, formal absolute inverse couplings",
              count, mass2.str(15), *[x.str(10) for x in absolute])
    print("One fully activated Dirac 15, supremum of projected correction", ceiling.str(16))


if __name__ == "__main__":
    verify()
