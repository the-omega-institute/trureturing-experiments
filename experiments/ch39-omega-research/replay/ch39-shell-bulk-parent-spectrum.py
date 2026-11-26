"""Exact indices, branching and face ordering for an SU(53) x SU(41) bulk spectrum."""
from collections import Counter
from fractions import Fraction


# Entries record complex dimension, circle charge and defining-normalized indices.
blocks = {
    "E": (12, 1, 3, 4, 0, 0),
    "J": (1, 6, 0, 0, 0, 0),
    "J*": (1, -6, 0, 0, 0, 0),
    "adjA": (15, 0, 8, 0, 0, 0),
    "adjB": (8, 0, 0, 6, 0, 0),
    "adjL": (3, 0, 0, 0, 4, 0),
    "adjR": (3, 0, 0, 0, 0, 4),
    "L": (2, 0, 0, 0, 1, 0),
    "R": (2, 0, 0, 0, 0, 1),
}
shared = Counter(adjA=1, adjB=1, adjL=2, adjR=2, L=1, R=1)
large = shared + Counter({"E": 1, "J*": 2})
small = shared + Counter({"J": 1, "J*": 1})


def invariants(branch):
    dim = det = center = 0
    indices = [0, 0, 0, 0]
    for name, count in branch.items():
        d, q, *quadratic = blocks[name]
        dim += count * d
        det += count * d * q
        center += count * d * q*q
        for i, value in enumerate(quadratic):
            indices[i] += count * value
    return dim, det, Fraction(center, 144), indices


d1, det1, c1, a1 = invariants(large)
d2, det2, c2, a2 = invariants(small)
assert (d1, d2, det1, det2) == (53, 41, 0, 0)
assert (c1, c2) == (Fraction(7, 12), Fraction(1, 2))
assert a1 == [11, 10, 9, 9] and a2 == [8, 6, 9, 9]
metric = [a + b + direct for a, b, direct in zip(a1, a2, (0, 0, 1, 1))]
assert metric == [19, 16, 19, 19]
assert [a-b for a, b in zip(a1, a2)] == [3, 4, 0, 0]
assert c1-c2 == Fraction(1, 12)
signed = large.copy()
signed.subtract(small)
assert {name: n for name, n in signed.items() if n} == {"E": 1, "J*": 1, "J": -1}
assert Fraction(d1-d2, 12) == 1


def traceless(n, j):
    rem = j % n
    return Fraction(rem*(n-rem), n)


def energies(cap, w4, w3, central):
    a, b = (1, 0) if cap else (-1, -1)
    return {k: w4*traceless(4, k-a) + w3*traceless(3, -k-b) + central*k*k
            for k in range(-5, 6, 2)}


# Both ends of the cone of positive high-energy SU couplings select the same faces.
for w4, w3, central in [(a1[0], a1[1], c1), (a2[0], a2[1], c2),
                        (metric[0], metric[1], c1+c2)]:
    for cap, optimum in ((True, -3), (False, -1)):
        values = energies(cap, w4, w3, central)
        if (w4, w3) == (8, 6) and cap:
            # At the isolated small-factor ray c=w3/12. Any large-factor term lifts this tie.
            assert sorted(k for k,v in values.items() if v == min(values.values())) == [-3, 1]
        else:
            assert [k for k,v in values.items() if v == min(values.values())] == [optimum]
assert Fraction(3, 8) == 1 / (1 + 1 + Fraction(2, 3))
print("bulk degrees, determinant charges:", (d1,d2), (det1,det2))
print("equal weights: SU4, family SU3, left, right levels", metric, "central", c1+c2)
print("target net index per SU4 component", -(2*(-3)+6*(-1))//4)
print("remaining complex spectator on the left before mass pairing: (J* - J) times weak L")
