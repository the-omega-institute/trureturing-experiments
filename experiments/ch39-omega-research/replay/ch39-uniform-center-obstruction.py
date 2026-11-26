"""Whole-window center obstruction with independent Wilson masses in both representations."""
from pathlib import Path
from runpy import run_path
from fractions import Fraction
from flint import arb, arb_poly, ctx

base = run_path(str(Path(__file__).with_name("ch39-rough-overlap-index.py")))
# Each mass interval has integer endpoints in units of 1/10000.
fields = {
    "A": ([(3,1,0,4,2),(0,2,3,4,1),(0,1,2,3,4),(0,2,1,4,3)],
          [(2,1,0,1,2),(2,1,2,1,0),(2,3,0,0,3),(2,1,2,1,2)],
          [(34540,35000,0)], [(34540,34740,-1),(34745,34818,-2),(34823,35000,-3)]),
    "B": ([(1,3,4,2,0),(4,2,3,1,0),(4,3,1,0,2),(2,4,0,3,1)],
          [(0,1,1,2,0),(1,3,1,3,2),(1,3,3,0,1),(0,0,0,3,1)],
          [(34540,34700,1)], [(34740,34745,10),(34818,34823,11)]),
    "C": ([(0,2,1,4,3),(2,0,1,4,3),(3,1,0,2,4),(0,3,2,1,4)],
          [(2,3,1,1,1),(1,2,1,1,1),(3,0,3,1,1),(0,3,0,1,2)],
          [(34700,35000,5)], [(34740,34823,3)]),
}


def exclude_zero(poly, lo, hi, depth=0):
    # The variable is the mass displacement from 139/40.
    if abs(poly(lo.union(hi))) > arb(2)**-64:
        return
    assert depth < 20
    mid = (lo+hi)/2
    exclude_zero(poly, lo, mid, depth+1)
    exclude_zero(poly, mid, hi, depth+1)


def verify_rep(values, intervals):
    d = values[0].nrows()
    chirality = base["stencil"](values)[0]*(arb(40)/61)
    windows = [(arb(low-34750)/10000, arb(high-34750)/10000) for low, high, _ in intervals]
    totals = [0]*len(intervals)
    for h in base["fibers"](values):
        # X(m)=gamma5 H(m0) - (m-m0) I. Gamma5-Hermiticity makes the polynomial real.
        determinant = arb_poly([c.real for c in (chirality*h).charpoly().coeffs()])
        for j, (lo, hi) in enumerate(windows):
            exclude_zero(determinant, lo, hi)
            center = h-chirality*((lo+hi)/2)
            polynomial = arb_poly([c.real for c in center.charpoly().coeffs()])
            cut = arb(2)**-35
            npos = base["positive_roots"](polynomial, cut)
            assert npos == base["positive_roots"](polynomial, -cut)
            totals[j] += 2*d-npos
    assert totals == [index for _, _, index in intervals]
    return totals


def check_cover():
    edges = sorted({x for _, _, bb, tt in fields.values() for a, b, _ in bb+tt for x in (a, b)})
    points = sorted(set(edges) | {Fraction(a+b, 2) for a, b in zip(edges, edges[1:])})
    for x in points:
        for y in points:
            assert any(a <= x <= b and c <= y <= d and (-i+2*j) % 5
                       for _, _, bb, tt in fields.values()
                       for a, b, i in bb for c, d, j in tt)
    assert arb(34540)/10000 < (25-5*arb(5).sqrt())/4
    assert 7**39 < 2**110


def verify():
    ctx.prec = 320
    check_cover()
    for name, (rows, exponents, bb, tt) in fields.items():
        values = base["field"](rows, exponents)
        left = verify_rep([u.conjugate() for u in values], bb)
        right = verify_rep([base["wedge"](u) for u in values], tt)
        print("Field, indices on certified intervals:", name, left, right, flush=True)


if __name__ == "__main__":
    verify()
