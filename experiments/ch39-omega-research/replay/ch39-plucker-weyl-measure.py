"""Flipped hypercharge, Plucker phase winding and the two-loop cost of the spectator U(1)."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from runpy import run_path
import sympy as sp
from flint import arb

source = run_path(str(Path(__file__).with_name("ch39-monopole-free-matching.py")))
model, q = source["model"], source["q"]


def classify():
    x = sp.symbols("x")
    a = 1/(1+6*x)
    up, down = sp.Rational(-2, 3), sp.Rational(1, 3)
    from_ten = a*(x+up)
    candidates = set(sp.solve((from_ten-up)*(from_ten-down), x))
    assert candidates == {-1, 0}
    for charge in candidates:
        scale = a.subs(x, charge)
        tens, fives, ones = charge, -3*charge, 5*charge
        physical = [scale*(v+tens) for v in (sp.Rational(1, 6), up, 1)]
        physical += [scale*(v+fives) for v in (down, sp.Rational(-1, 2))]
        physical += [scale*ones]
        expected = ([sp.Rational(1, 6), down, 0, up, sp.Rational(-1, 2), 1]
                    if charge == -1 else [sp.Rational(1, 6), up, 1, down, sp.Rational(-1, 2), 0])
        assert physical == expected
        assert 3*tens+fives == 10*tens+5*fives+ones == 10*tens**3+5*fives**3+ones**3 == 0
    assert [(q5+2*qx) % 5 for q5, qx in ((2,-1),(-1,3),(0,-5))] == [0]*3
    for i5, i10 in product(range(5), repeat=2):
        center = -i5+2*i10
        assert (3*i5-i10-2*center) % 5 == 0
        if center % 5:
            assert (3*center) % 5 != 0
    assert sum(dim*charge**2 for dim, charge in [(10,-1),(5,3),(1,-5)]) == 80


def kinetic(gauge):
    pairs = []
    def observe(j, inv, deriv):
        if j >= 23*64 and 0 in (inv[1]-q("16/15")*inv[2]):
            g_x = q("5/3")*(25*inv[0]-q("24/25")*inv[2])
            pairs.append((g_x, g_x/inv[2]))
    model["running"]["integrate"](gauge, 33*64, observe)
    assert pairs
    return pairs


def verify():
    classify()
    centers = model["GAUGE"]
    for g_x, ratio in kinetic([q(v) for v in centers]):
        assert q("126.94") < g_x < q("127.12") and q("38.53") < ratio < q("38.67")
    for sides in product((0, 1), repeat=3):
        gauge = [q(c*(1+F(k-1, 1000))).union(q(c*(1+F(k, 1000)))) for k, c in zip(sides, centers)]
        for g_x, ratio in kinetic(gauge):
            assert q("126.27") < g_x < q("127.79") and q("38.20") < ratio < q("39.00")
            assert 1/(4*arb.pi()*g_x) < q("0.000631")
    print("Two hypercharge assignments. Flipped X kinetic coefficient 126.27 to 127.79.")
    print("The diagonal center is neutral in the flipped assignment; Plucker disk winding persists.")


if __name__ == "__main__":
    verify()
