"""Exact face counts and hinge orders for finite four-dimensional Regge comparisons."""
from collections import Counter
from itertools import combinations

F = frozenset


def one_five(complex_, facet, center):
    assert facet in complex_ and all(center not in v for v in complex_)
    replacements = {F((center,)) | (facet-{removed}) for removed in facet}
    return (complex_-{facet}) | replacements


def two_four(complex_, face):
    incident = {v for v in complex_ if face <= v}
    assert len(face) == 4 and len(incident) == 2
    nodes = F.union(*incident)
    opposite = nodes-face
    assert len(opposite) == 2 and not any(opposite <= v for v in complex_)
    return (complex_-incident) | {nodes-{removed} for removed in face}


def hinges(complex_):
    return Counter(F(face) for simplex in complex_ for face in combinations(sorted(simplex),3))


def counts(complex_):
    tetra = Counter(F(face) for simplex in complex_ for face in combinations(sorted(simplex),4))
    assert all(v==2 for v in tetra.values())
    orders = hinges(complex_)
    assert sum(orders.values())==10*len(complex_)
    return len(complex_), len(orders), sum(q*q for q in orders.values())


start = {F(v) for v in combinations(range(6),5)}
first = one_five(start, F(range(5)), 6)
second = two_four(first, F(range(4)))
rows = [counts(T) for T in (start, first, second)]
assert [r[:2] for r in rows] == [(6,20),(10,30),(12,34)]
# Break the facet symmetry by one specified further subdivision.
third = one_five(second, min(second,key=lambda v: tuple(sorted(v))), 7)
alternatives = [(counts(one_five(third,facet,8)), tuple(sorted(facet))) for facet in third]
low, high = min(alternatives), max(alternatives)
assert low[0][:2] == high[0][:2] == (20,54) and low[0][2] < high[0][2]
print("N4,N2,sum q^2",rows,counts(third))
print("alternative one-five moves",low,high)
print("difference in sum q^2",high[0][2]-low[0][2])
