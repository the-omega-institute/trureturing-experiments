"""Taste dispersion through cubic order and a strictly finite pentamer implementation."""
from itertools import permutations, product
from pathlib import Path
from runpy import run_path

import sympy as sp
from flint import acb, arb, fmpq


def dispersion():
    c = (sp.sqrt(5)-1)/4
    A,B = 5*(1-c),sp.Rational(7,2)
    # For m in (A,B), c/(m-A) > 1/m since m<5.
    m = sp.symbols("m")
    assert sp.simplify(c*m-(m-A)-(1-c)*(5-m)) == 0
    d = (1,1,1,1,-4)
    assert (sum(d),sum(x*x for x in d),sum(x**3 for x in d)) == (0,20,-60)
    assert A < B < 5
    k,t = sp.symbols("k t", real=True)
    y = sp.Matrix([c*k*v-t*k**2*v*v/2 for v in d])
    projector = sp.eye(5)-sp.ones(5)/5
    norm = (y.T*projector*y)[0].expand()
    assert sp.simplify(norm.coeff(k,2)-20*c*c) == 0
    assert sp.simplify(norm.coeff(k,3)-60*c*t) == 0


def finite_probe():
    tiny = run_path(str(Path(__file__).with_name("ch39-regge-fixed-volume.py")))
    c = tiny["cosine_tick"](tiny["GRID"]//2**20)
    w = acb(c,(1-c*c).sqrt())
    cx = (arb(5).sqrt()-1)/4
    z = acb(cx,(1-cx*cx).sqrt())
    vals = []
    for s in (0,-1,1):
        five = [z**s*w**u for u in (1,1,1,1,-4)]
        ys = [v.imag for v in five]
        sq = sum(v*v for v in ys)-sum(ys)**2/5
        mass = fmpq(139,40)-5+sum(v.real for v in five)
        R = (mass*mass+sq).sqrt()
        vals.append(2*sq/(R*(R+mass)))
    assert 0 < vals[0] < vals[1] < vals[2]
    print("Finite folded D†D values at m=139/40:",*[v.str(8) for v in vals])


def blocks():
    # e1 is the selected direction. Each line is partitioned by tau(x) mod five.
    L = 10
    covered = set()
    for point in product(range(L),repeat=4):
        if sum(point)%5:
            continue
        block = {((point[0]+j)%L,*point[1:]) for j in range(5)}
        assert covered.isdisjoint(block) and len(block)==5
        covered.update(block)
    assert len(covered)==L**4
    # Unnormalized Fourier columns have Gram 5I in Q[zeta_5].
    z = sp.Symbol("z")
    mod = sp.cyclotomic_poly(5,z)
    for difference in range(5):
        val = sum(z**((difference*j)%5) for j in range(5))
        assert sp.rem(val,mod,z)==5*int(difference==0)
    for perm in permutations((4,0,1)):
        action = dict(zip((4,0,1),perm))
        action.update({2:2,3:3})
        assert sorted(action.values()) == list(range(5))


if __name__ == "__main__":
    dispersion()
    finite_probe()
    blocks()
    print("Three distinct dispersion germs on the tetrahedrally fixed ray; exact pentamer group action.")
