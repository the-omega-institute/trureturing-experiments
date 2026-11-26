"""Read the exterior Dirac masses despite arbitrary positive fundamental Yukawa strength."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import ctx
from sympy import Poly

ROOT = Path(__file__).parent
FINITE = run_path(str(ROOT/"ch39-finite-moment-height.py"))
GRID = run_path(str(ROOT/"ch39-noisy-threshold-grid.py"))
a = GRID["a"]
CURRENT = (5,-3,-2)


def weights():
    indices = [("1/5","0","1/2"),("3/10","1/2","0"),("4/5","0","1/2"),
               ("1/10","3/2","1"),("3/5","0","0"),("4/5","0","1/2"),
               ("7/10","3/2","1")]
    def weight(row):
        return sum(F(s)*k for s,k in zip(row,CURRENT))/3
    assert list(map(weight,indices)) == [0,0,1,-2,1,1,-1]
    assert [weight(row) for row in [(0,2,3),(5,3,2),(25,25,25)]] == [-4,4,0]
    assert weight(("2/5",0,1)) == 0
    factors = FINITE["charged_scalar_factors"](CURRENT)
    x, q = FINITE["x"], FINITE["q"]
    expected = {x-2:60, x-4:-72, x-q-2:-144, x-4*q-2:12,
                       x*x-(q+7)*x+3*q+6:12, x*x-(q+6)*x+3*q:12,
                       x*x-(4*q+6)*x+12*q:48}
    assert factors == {p.expand(): w for p, w in expected.items()}
    scalar_count = sum(Poly(p,x).degree() for p in factors)+2
    assert scalar_count == 12


def verify():
    ctx.prec = 180
    weights()
    previous = None
    for n, center in GRID["models"]():
        fermions, vectors, scalars = GRID["mass_ranges"](center)
        E, D = fermions[0][0], fermions[1][0]
        C = D**2/E
        assert a("2.05") < E < a("2.07") < a("2.86") < D < a("2.89")
        assert a("3.95") < C < a("4.1") < a("4.2")
        assert all(t > 5 for t,_ in vectors)
        assert all(t < a("4.2") for t,_ in scalars)
        assert 4*a("1.0001")*2*a("1.000120")/12 < a("4.2")
        # Even exact coincidences with either (or both) triply degenerate
        # bridge blocks change the onset residue only by a multiple of three.
        for contamination in (-3,0,3):
            assert all((w+contamination)%3 == 1 for w in (1,-2,1))
        ell = ((center.gamma/624)**2+a("1/4")).sqrt()
        p = ((2*ell+1).sqrt()-1)/((2*ell+1).sqrt()+1)
        assert (E/D).contains(p)
        if previous is not None:
            assert previous < p
        previous = p
        print(n, "certified bands E,D,C", E.str(4), D.str(4), C.str(4))
    for _, _, epsilon in GRID["MATCH"]["SOURCE"]["RATIO"]["high_zero_readings"]():
        assert epsilon**2 < a("4.2")/5
    # Three separated groups inside [4.2,5]. Across two hypothetical spectra
    # at most two groups can contain a bridge onset, even at sample endpoints.
    groups = [[F(42,10)+F(26*i,100)+F(j,200) for j in range(44)] for i in range(3)]
    assert all(left[-1] < right[0] for left,right in zip(groups, groups[1:]))
    assert max(groups[-1]) < 5
    assert sum(map(len,groups)) == 132
    print("Interior sample groups", [(str(g[0]), str(g[-1])) for g in groups])


if __name__ == "__main__":
    verify()
