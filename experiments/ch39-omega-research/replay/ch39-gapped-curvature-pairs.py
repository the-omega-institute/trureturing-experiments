"""Paired singlet and triplet matter at a common mass offset. All scales use Q or M=Q/r."""
from fractions import Fraction
from pathlib import Path
import runpy
from flint import arb

root = Path(__file__).resolve().parent
graph = runpy.run_path(str(root / "ch39-neutral-clock-consecutive.py"))
q, r = graph["q"], graph["r"]
pi = arb.pi()
rotations, _ = graph["rotations"](719)
capacity = graph["pot"]["short"]["base"]["arc_budget"]()[1].upper()
degree = capacity / 2

def shifted_moment(mass_squared):
    # Shift is in units of M^2. Positive walk weights decrease with Z.
    lo = mass_squared
    lower = q(0)
    upper = q(0)
    step = q(2)**(q(1)/16)
    for _ in range(32 * 16):
        hi = lo * step
        lower += (hi - lo) * graph["point"](hi, False, rotations).lower()
        upper += (hi - lo) * graph["point"](lo, True, rotations).upper()
        lo = hi
    assert lo > 2 * degree
    upper += 608 * degree**3 / lo**2
    return lower.lower().union(upper.upper()) / r**2

def check():
    # Eight unshifted real scalars and ten unshifted two-component fermions
    # accompany the vector singlet. The second vector heat is the triplet.
    assert 8 - 2 * 10 + 3 == -9
    assert 8 * (Fraction(1, 6) - Fraction(1, 8)) + Fraction(10, 6) - Fraction(1, 2) == Fraction(3, 2)
    assert 2 * (Fraction(1, 6) - Fraction(1, 4)) + Fraction(1, 6) == 0
    assert 10 + 3 + 3 == 16
    assert 8 + 6 + 3 * 2 == 20
    shifted = shifted_moment(q(1))
    unshifted = q(".000048797").union(q(".00030629"))
    charged = q(".05").union(q("1.99779"))
    scale = (1000 * q("26.8220").exp()).union(1000 * q("27.0455").exp())
    planck = q("1.220890e19")
    # Interval division with independent positive endpoints avoids dependency
    # overestimates from reciprocal at an interval spanning several decades.
    need = pi * (planck / scale)**2 - 2 * charged / 3 - unshifted / 4
    a = (need.lower() / shifted.upper()).union(need.upper() / shifted.lower())
    weyl = runpy.run_path(str(root / "ch39-curvature-spectrum-rigidity.py"))["weyl"]
    # Six shifted real directions at each coupling, eight at 1/8, four neutral
    # vector directions (Maxwell plus compensating scalar at its zero mode).
    # The charged sector has rank two relative to one ordinary copy.
    neutral_const = q(6)/144 + q(4)/576 + q(4)/72
    charged_const = 2 * (q(5)/900 + q(2)/72)
    constant = neutral_const + charged_const
    slope_lo = weyl * (6*a.lower()**2 + constant)/(16*pi**2)
    slope_hi = weyl * (6*a.upper()**2 + constant)/(16*pi**2)
    # Ordinary shifted scalar has an internal constant mode. On a four-sphere
    # the external constant mode hits zero when R = M^2/(a-1/4).
    rc = ((scale.lower()/r)**2/(a.upper()-q(1)/4)).union(
         (scale.upper()/r)**2/(a.lower()-q(1)/4))
    ir = 4*Fraction(1,24)**2 + Fraction(1,72)
    assert ir == Fraction(1,48)
    assert q(".0000367549") < shifted < q(".0000421204")
    assert q("3.4733e13") < a < q("6.4254e13")
    assert q("3.5959e27") < slope_lo < slope_hi < q("1.2306e28")
    assert q("1.3888e5") < rc.sqrt() < q("2.3620e5")
    for label, val in (("damped moment", shifted), ("fitted a", a),
                       ("R^2 UV slope / (Lambda/Q)", slope_lo.lower().union(slope_hi.upper())),
                       ("sqrt(Rcrit) GeV", rc.sqrt())):
        print(label, val.lower().str(13), val.upper().str(13))
    print("normalized IR heat coefficient of pure R^2", ir)

if __name__ == "__main__":
    check()
