"""Uniform derivative certificate for the original cubic shell's top-three sum."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx
from sympy import Poly, Rational, prod, symbols

shell = run_path(str(Path(__file__).with_name("ch39-refined-mass-obstruction.py")))
ball = shell["ball"]


def verify():
    ctx.prec = 224
    nodes, edges = shell["geometry"]()
    charged = [F(0)]*len(nodes)
    for i,j,a,k in edges:
        if k:
            charged[i] += a
            charged[j] += a
    # Both row and column sums of |M(t)-M(0)| for complex |t| <= 1/2.
    delta = 2*max(charged)
    assert delta < F(21,10000)
    sample = [F(j,112)**2 for j in range(5)]
    readings = []
    for j in range(5):
        it = acb(0, ball(F(j,112)))
        phase = (1+it)/(1-it)
        matrix = shell["matrices"](nodes, edges, phase)
        eig = matrix.eig(multiple=True)
        upper = [v.real for v in eig if v.real > arb(".075")]
        lower = [v.real for v in eig if v.real < arb(".075")]
        assert len(upper) == 3 and len(lower) == len(nodes)-3
        if j == 0:
            assert min(upper) > arb(".0795") and max(upper) < arb(".0796")
            assert max(lower) < arb(".070")
        readings.append(sum(upper))
    # The perturbed upper cluster on the complex disk lies within delta of the
    # unperturbed three levels; its sum has modulus less than 3*(.0796+delta).
    bound = F(1,4)
    assert 3*(F(".0796")+delta) < bound and 2*delta < F(".0795")-F(".070")
    s = symbols("s")
    coeffs = [arb(0) for _ in range(4)]
    for j, reading in enumerate(readings):
        basis = prod((s-Rational(v))/Rational(sample[j]-v) for i,v in enumerate(sample) if i != j)
        poly = Poly(basis.diff(s), s)
        for k in range(4):
            coeffs[k] += reading*ball(F(poly.nth(k)))
    interval = arb(0).union(ball(sample[-1]))
    derivative = arb(0)
    for c in reversed(coeffs):
        derivative = derivative*interval+c
    radius, width = F(1,4), sample[-1]
    dist = radius-width
    remainder = radius*bound/dist**5*(5*width**4/dist + width**5/dist**2)
    assert remainder < F("3.7e-9")
    enclosure = derivative + arb(0, ball(remainder))
    assert arb("-7.5e-7") < enclosure < arb("-7.2e-7")
    assert acb.zeta_zero(1).imag > 14 and width == F(1,28)**2
    print("Lagrange derivative",derivative.str(12),"uniform remainder",float(remainder))
    print("Full top-three derivative on [0,1/784]",enclosure.str(12))


if __name__ == "__main__":
    verify()
