"""Index the low heavy-scalar spectrum on the graph carrying the flavor clocks."""
from fractions import Fraction
from functools import cache
from pathlib import Path
import runpy
from types import FunctionType

from flint import acb, acb_mat, arb

root = Path(__file__).resolve().parent
levels = runpy.run_path(str(root / "ch39-quaternion-scalar-levels.py"))
shell = levels["src"]
height = acb.zeta_zero(719).imag
p, q = [shell[key] * height for key in ("a", "b")]
cos_axis = levels["c"]
sin_axis = levels["v"]
radius = arb("2e-7")

def clock_count(width, all_sectors=False):
    su, cu = (p + arb(0, width)).sin_cos()
    sv, cv = (q + arb(0, width)).sin_cos()
    U = acb_mat([[acb(cu, su), 0], [0, acb(cu, -su)]])
    V = acb_mat([[acb(cv, cos_axis * sv), sin_axis * sv],
                 [-sin_axis * sv, acb(cv, -cos_axis * sv)]])
    fn = levels["count"]
    context = dict(fn.__globals__, U=U, V=V)
    if all_sectors:
        context["sectors"] = [acb_mat([[int(i == j) for j in range(38)] for i in range(38)])]
    return FunctionType(fn.__code__, context)

count = clock_count(0)

@cache
def index_at(energy):
    return count(levels["q"](Fraction(energy)))

def isolate():
    assert index_at(Fraction(1, 10)) == 0 and index_at(16) == 18
    assert index_at(446**2) == 18
    rows = []
    for target in range(1, 19):
        lo, hi = Fraction(1, 10), Fraction(16)
        for _ in range(21):
            mid = (lo + hi) / 2
            if index_at(mid) < target:
                lo = mid
            else:
                hi = mid
        rows.append((lo, hi))
    return rows

def window_levels():
    rows = isolate()
    full = clock_count(radius, all_sectors=True)
    intervals = []
    for pair in range(9):
        lo, hi = rows[2*pair]
        assert rows[2*pair + 1] == (lo, hi)
        # Recheck both ends throughout the clock square without cyclic reduction.
        a = levels["q"](lo) - arb(".001")
        b = levels["q"](hi) + arb(".001")
        assert full(a) == 2*pair and full(b) == 2*pair + 2
        intervals.extend([a.union(b)] * 2)
    assert full(arb(446**2)) == 18
    return intervals

if __name__ == "__main__":
    intervals = window_levels()
    for i, (lo, hi) in enumerate(isolate()[::2]):
        print(i+1, float(lo), float(hi), intervals[2*i].str(10))
