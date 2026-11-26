"""Extrapolate the two fitted shell clocks to both fattest six-digit cube shells."""
from collections import defaultdict
from fractions import Fraction
from itertools import permutations, product
from math import comb
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

ctx.prec = 96
penalty = arb(1) / 51**5
choices = ((26,), (29,), (26, 29))


def dot(v, w):
    return sum(x*y for x, y in zip(v, w))


def cross(v, w):
    return tuple(v[(i+1)%3]*w[(i+2)%3] - v[(i+2)%3]*w[(i+1)%3] for i in range(3))


def fattest(m):
    levels = defaultdict(list)
    for v in product(range(m), repeat=3):
        levels[dot(v, v)].append(v)
    size = max(map(len, levels.values()))
    return {r: levels[r] for r in sorted(levels) if len(levels[r]) == size}


def clock(v, w):
    r = dot(v, v)
    assert r == dot(w, w)
    edges = sorted((r-dot(v, p), p) for p in set(permutations(w)))
    lengths = [d for d, p in edges]
    assert len(set(lengths)) == len(lengths)
    d0, d1, d2 = lengths[:3]
    normals = [cross(v, p) for _, p in edges[1:3]]
    u, z = normals
    cosine2 = Fraction(dot(u, z)**2, dot(u, u)*dot(z, z))
    assert dot(u, z) > 0
    cosine = (arb(cosine2.numerator)/cosine2.denominator).sqrt()
    angles = [(1-arb(d)/r).acos() for d in (d0, d1, d2)]
    return (d0, d1, d2), cosine, (angles[1]-angles[0], angles[2]-angles[0])


base = {17: clock((0, 1, 4), (2, 2, 3)), 867: clock((1, 5, 29), (13, 13, 23))}
extended = {}
for r, points in fattest(6).items():
    orbits = sorted({tuple(sorted(v)) for v in points})
    assert len(orbits) == 2 and orbits[0][0] == 0 < orbits[1][0]
    extended[r] = clock(*orbits)
shells = base | extended
assert [shells[r][0] for r in (17, 867, 26, 29)] == [(3,6,7),(122,362,402),(3,5,7),(3,5,6)]


def readout(c, p, q):
    s, x = p.sin_cos()
    t, y = q.sin_cos()
    return (s*s, t*t, 1-(x*y+c*s*t)**2)


def measured(h, r):
    _, c, (f, g) = shells[r]
    return readout(c, h*f, h*g)


def reference(r, i, j):
    return readout(shells[r][1], arb.pi()*i/9, arb.pi()*j/9)


def distance(x, y):
    return sum((a-b)**2 for a,b in zip(x,y))


original = {17: reference(17, 2, 1), 867: reference(867, 4, 6)}
target = {(r, a): reference(r, a, 1 if r == 26 else 2) for r in extended for a in range(9)}


def pieces(h):
    old = sum(distance(measured(h, r), ref) for r, ref in original.items())
    new = {r: measured(h, r) for r in extended}
    extra = {(r,a): distance(new[r], ref) for (r,a),ref in target.items()}
    return old + penalty*h*h, extra


def total(old, extra, key):
    layers, a = key
    return old + sum(extra[r,a] for r in layers)


def main():
    for m in range(3,7):
        levels = fattest(m)
        print("digits", m, "largest size", len(next(iter(levels.values()))), "squared radii", list(levels))
    quartic = [(3,-4), (4,-2)]
    local = run_path(str(Path(__file__).with_name("ch39-shell-local-nine-flux.py")))
    witnesses = [quartic] + [row[1] for row in local["solve"](54)]
    def labels(coeff):
        return tuple(sum(b*(comb(d,k)-comb(3,k)) for k,b in coeff)%9 for d in (5,6,7))
    assert list(map(labels, witnesses)) == [(8,2,1),(5,2,1),(5,2,1)]
    print("new and inherited relative labels", list(map(labels, witnesses)))

    # Seed upper bounds and certify uniqueness with ordered Arb intervals at each zero.
    winners = {}
    energies719 = {}
    keys = [(layers, a) for layers in choices for a in range(9)]
    def compare(n, h, keys):
        old, extra = pieces(h)
        for key in keys:
            value = total(old, extra, key)
            if n == 719:
                energies719[key] = value
            if key not in winners or value < winners[key][1]:
                winners[key] = n, value
            else:
                assert value > winners[key][1]
    for n,z in enumerate(acb.zeta_zeros(1,1024), start=1):
        compare(n, z.imag, keys)
    start = int(acb.zeta_zero(1024).imag.floor().unique_fmpz())
    cap = max(int((v/penalty).sqrt().ceil().unique_fmpz()) for n,v in winners.values())
    assert all(penalty*cap**2 > v for n,v in winners.values())

    def cover(lo, hi, keys):
        h = lo.union(hi)
        old, extra = pieces(h)
        possible = [key for key in keys if not total(old, extra, key) > winners[key][1]]
        if not possible:
            return
        if hi-lo > arb(1)/8:
            mid = (lo+hi)/2
            cover(lo, mid, possible)
            cover(mid, hi, possible)
            return
        nl, nh = (v.zeta_nzeros().unique_fmpz() for v in (lo,hi))
        assert nl is not None and nh is not None
        for n in range(max(1024,int(nl))+1, int(nh)+1):
            compare(n, acb.zeta_zero(n).imag, possible)
    for h in range(start,cap):
        cover(arb(h), arb(h+1), keys)
    expected = ([2619,2619,719,719,651,651,651,719,2619],
                [583,2447,2454,2619,719,719,719,858,357],
                [4372,4372,2454,6621,4867,2166,2121,719,4372])
    for layers, want in zip(choices, expected):
        assert [winners[layers,a][0] for a in range(9)] == want
        print("shells", layers, "all nine unique winners", [winners[layers,a][0] for a in range(9)])
        for a in (8,5):
            n,v = winners[layers,a]
            print("phase", a, "winner", n, "minimum", v.str(14),
                  "719 excess", (energies719[layers,a]-v).str(14))
    print("uniform infinite tail starts at height", cap)


if __name__ == "__main__":
    main()
