"""Exact rational mass response at the first cyclic cutoff of the shell graph."""
from collections import Counter
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx
from sympy import Matrix, Poly, Rational, cancel, eye, fraction, symbols

U = symbols("u")
T = Rational(1, 19)


def response(power):
    a, b, c = [Rational(12, d)**power for d in (32, 52, 56)]
    p2 = 1+2*b*b+2*c*c
    companion = Matrix([[0, 0, 0, -4*b*b*c*c*(1-U*U)],
                        [1, 0, 0, 4*b*c*U], [0, 1, 0, p2], [0, 0, 1, 0]])
    block = T*companion
    den = eye(4)-block**3
    return cancel(3*(block*den.adjugate()).trace()/den.det()+6*(T*a)**4/(1-(T*a)**6))


def positive_on_unit_interval(poly):
    p = Poly(poly, U)
    assert p.eval(0) > 0 and p.eval(1) > 0
    assert p.count_roots(0, 1) == 0


def evaluated(poly, arg):
    val = arb(0)
    for coeff in Poly(poly, U).all_coeffs():
        val = val*arg+arb(int(coeff.p))/int(coeff.q)
    return val


def graph_blocks():
    points = [p for p in product(range(30), repeat=3) if sum(x*x for x in p) == 867]
    blocks = [{i} for i in range(len(points))]
    edges = {}
    for (i, p), (j, q) in combinations(enumerate(points), 2):
        d = sum((a-b)**2 for a, b in zip(p, q))
        if d <= 56:
            edges[i, j] = d
            parts = [component for component in blocks if {i,j} & component]
            blocks = [component for component in blocks if component not in parts] + [set.union(*parts)]
    assert Counter(map(len, blocks)) == {1: 1, 2: 3, 4: 3}
    for component in blocks:
        inside = {(i,j): d for (i,j), d in edges.items() if i in component}
        if len(component) == 4:
            assert Counter(inside.values()) == {12: 1, 52: 2, 56: 2}
            ends, = [ij for ij, d in inside.items() if d == 12]
            types = [sorted(d for ij, d in inside.items() if i in ij) for i in ends]
            assert sorted(types) == [[12, 52, 52], [12, 56, 56]]
        elif len(component) == 2:
            assert list(inside.values()) == [32]


def verify():
    ctx.prec = 180
    graph_blocks()
    model = run_path(str(Path(__file__).with_name("ch39-shell-zero-mass-hierarchy.py")))
    a, _, fit = model["parameters"]()
    for p in range(3):
        num, den = fraction(cancel(response(p)-Rational(46913, 10**9)))
        positive_on_unit_interval(num*den)
    num, den = fraction(response(3))
    positive_on_unit_interval(num.diff(U)*den-num*den.diff(U))
    e0, e1 = (evaluated(num, arb(u))/evaluated(den, arb(u)) for u in (0, 1))
    low, v, log_x, g5, g10 = fit(e0.union(e1), a(1))
    assert v*log_x.exp()*e0 > 1000
    assert a(".618") < low < a(".621") and g5**2 < arb.pi() and a("1.33") < log_x.exp() < a("1.35")
    assert a("7.59e14") < v < a("7.63e14")
    readings = []
    for j in range(11485, 11491):
        zero = acb.zeta_zero(j)
        assert zero.real == a("1/2")
        u = 1-1/(2*(zero.imag**2+a("1/4")))
        readings.append(evaluated(num, u)/evaluated(den, u))
    gaps = [y-x for x,y in zip(readings, readings[1:])]
    assert all(a("7e-24") < gap < a("1.8e-23") for gap in gaps)
    print("Cubic weight: g_low", low.str(8), "v/GeV", v.str(9), "x", log_x.exp().str(8))
    print("Five certified differences:", [x.str(10) for x in gaps])


if __name__ == "__main__":
    verify()
