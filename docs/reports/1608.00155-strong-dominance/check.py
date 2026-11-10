#!/usr/bin/env python3
"""Exact independent checks for arXiv:1608.00155v1's strong-dominance conjecture.

Only Python's standard library is needed. Every asserted value uses integers
or Fraction. Polynomial arithmetic is over Q; no floating point is used.
"""
from fractions import Fraction as Q
from itertools import combinations, product
from math import comb, prod


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for tail in compositions(total - first, length - 1):
                yield (first,) + tail


def word_counts(weights, n):
    """Numerators of the iid-word PMF; denominator sum(weights)**n."""
    out = [0] * (len(weights) + 1)
    for w in product(range(len(weights)), repeat=n):
        out[len(set(w))] += prod(weights[j] for j in w)
    assert sum(out) == sum(weights) ** n
    return out


def source_pmf(weights, n, k):
    """Source Section 3: inclusion-exclusion over all subsets A of H."""
    N = len(weights)
    return sum((-1) ** (k - j) * comb(N - j, k - j)
               * sum(sum(weights[i] for i in A) ** n
                     for A in combinations(range(N), j))
               for j in range(k + 1))


def majorized(p, q):
    a, b = sorted(p, reverse=True), sorted(q, reverse=True)
    return sum(a) == sum(b) and all(sum(a[:k]) <= sum(b[:k])
                                   for k in range(1, len(a)))


def single_crossing(delta):
    """A real threshold exists iff no positive entry precedes a negative one."""
    seen_positive = False
    for x in delta:
        if x > 0:
            seen_positive = True
        if x < 0 and seen_positive:
            return False
    return True


# Ascending coefficient tuples implement exact Q[t].
def trim(a):
    a = list(map(Q, a))
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return tuple(a)


def add(a, b):
    return trim([(a[i] if i < len(a) else 0) +
                 (b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b)))])


def scale(c, a):
    return trim([c * x for x in a])


def mul(a, b):
    c = [Q(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return trim(c)


def power(a, n):
    c = (Q(1),)
    for _ in range(n):
        c = mul(c, a)
    return c


def psum(polys):
    out = (Q(0),)
    for a in polys:
        out = add(out, a)
    return out


def evaluate(a, t):
    return sum(c * t ** i for i, c in enumerate(a))


def interval_bounds(a, upper):
    """Exact coefficientwise enclosure for all 0 <= t <= upper."""
    return (a[0] + sum(min(Q(0), a[i] * upper ** i)
                       for i in range(1, len(a))),
            a[0] + sum(max(Q(0), a[i] * upper ** i)
                       for i in range(1, len(a))))


def word_polynomials(weights, n):
    out = [(Q(0),)] * (len(weights) + 1)
    for w in product(range(len(weights)), repeat=n):
        term = (Q(1),)
        for j in w:
            term = mul(term, weights[j])
        out[len(set(w))] = add(out[len(set(w))], term)
    return out


def counting_polynomials(weights):
    """Independent singleton/pair/multiplicity formulas at N=4, n=5."""
    one = psum(power(w, 5) for w in weights)
    pair_sum = psum(power(add(weights[i], weights[j]), 5)
                    for i, j in combinations(range(4), 2))
    two = add(pair_sum, scale(-3, one))
    four = (Q(60),)
    for w in weights:
        four = mul(four, w)
    three = add((Q(1),), scale(-1, psum([one, two, four])))
    return [(Q(0),), one, two, three, four]


# Fidelity: reproduce two exact source identities on 600 finite input cases.
fidelity_cases = 0
for N in range(1, 5):
    for w in compositions(6, N):
        for n in range(1, 6):
            numerators = word_counts(w, n)
            assert all(numerators[k] == source_pmf(w, n, k)
                       for k in range(min(N, n) + 1))
            # Source Section 3: E(Y_n) = N - sum_i (1-p_i)**n.
            assert sum(k * numerators[k] for k in range(N + 1)) == (
                N * 6 ** n - sum((6 - x) ** n for x in w))
            fidelity_cases += 1
print('FIDELITY inclusion-exclusion and E(Y_n):', fidelity_cases, 'exact cases')

# Source's stated n=2,3 and n=4,N<=5 result on a specified rational grid.
grid = list(compositions(6, 4))
small_pairs = 0
for n in [2, 3, 4]:
    pmfs = {w: word_counts(w, n) for w in grid}
    for p in grid:
        for q in grid:
            if majorized(p, q):
                assert single_crossing([a-b for a, b in zip(pmfs[p], pmfs[q])])
                small_pairs += 1
print('FIDELITY strong dominance: N=4,n=2,3,4,denominator=6;',
      small_pairs, 'majorized ordered pairs checked')

# Concrete counterexample, from independent enumeration of 4**5=1024 words.
p, q, n = (48, 48, 2, 2), (50, 46, 3, 1), 5
D = 100 ** n
pp, qq = word_counts(p, n), word_counts(q, n)
assert all(x > 0 for x in p + q) and majorized(p, q)
assert pp == [0, 509608000, 7874904000, 1560192000, 55296000]
assert qq == [0, 518463220, 7873379340, 1566757440, 41400000]
delta = [a-b for a, b in zip(pp, qq)]
assert delta == [0, -8855220, 1524660, -6565440, 13896000]
assert delta[2] > 0 and delta[3] < 0 and not single_crossing(delta)
cdf = [sum(delta[:k+1]) for k in range(1, 5)]
assert all(x < 0 for x in cdf[:-1]) and cdf[-1] == 0
print('WITNESS N=4 n=5 words=1024 D=', D)
print('p prefixes:', [sum(p[:k]) for k in range(1, 5)])
print('q prefixes:', [sum(q[:k]) for k in range(1, 5)])
print('P_p numerators k=1..4:', pp[1:])
print('P_q numerators k=1..4:', qq[1:])
print('Delta numerators k=1..4:', delta[1:])
print('CDF difference numerators k=1..4:', cdf)
print('NO_THRESHOLD: Delta_2>0 forces c<2, whereas Delta_3<0 requires c>=3')

# The full real interval is certified by exact polynomial identities and bounds.
pfam = [(Q(1, 2), -2), (Q(1, 2), -2), (0, 2), (0, 2)]
qfam = [(Q(1, 2),), (Q(1, 2), -4), (0, 3), (0, 1)]
u = Q(1, 100)
for fam in [pfam, qfam]:
    assert psum(fam) == (Q(1),)
    for a in fam:
        lo, _ = interval_bounds(a, u)
        assert lo > 0 or (a[0] == 0 and len(a) == 2 and a[1] > 0)
    for a, b in zip(fam, fam[1:]):
        assert interval_bounds(add(a, scale(-1, b)), u)[0] >= 0
for k in range(1, 5):
    assert interval_bounds(add(psum(qfam[:k]), scale(-1, psum(pfam[:k]))), u)[0] >= 0
assert [evaluate(a, u) for a in pfam] == [Q(x, 100) for x in p]
assert [evaluate(a, u) for a in qfam] == [Q(x, 100) for x in q]
P, QP = word_polynomials(pfam, 5), word_polynomials(qfam, 5)
assert P == counting_polynomials(pfam) and QP == counting_polynomials(qfam)
dpolys = [add(a, scale(-1, b)) for a, b in zip(P, QP)]
assert dpolys[2] == (0, 0, 5, -360, 1270, -2340)
assert dpolys[3] == (0, 0, -10, 360, -1670, 1560)
print('FAMILY 0<t<=1/100: positive sorted normalized vectors, p_t majorized by q_t')
for k in range(1, 5):
    a = dpolys[k]
    assert a[0] == a[1] == 0
    lo, hi = interval_bounds(a[2:], u)
    assert hi < 0 if k in [1, 3] else lo > 0
    assert evaluate(a, u) == Q(delta[k], D)
    print('Delta_'+str(k)+' coefficients (t^0 upward):', ','.join(map(str, a)))
    print('Delta_'+str(k)+'/t^2 enclosure:', str(lo), '<= bracket <=', str(hi))
assert psum(dpolys) == (Q(0),)
print('FAMILY sign pattern (-,+,-,+) certified on the full real interval')
print('ALL_CHECKS_OK')
