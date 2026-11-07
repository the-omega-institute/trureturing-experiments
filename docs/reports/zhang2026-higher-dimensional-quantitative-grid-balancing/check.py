#!/usr/bin/env python3
"""Exact finite checks for Zhang arXiv:2606.24016v2, Section 6.

All evidence uses integers and fractions. No floating point is used.
A reduced combinatorial Laplacian counts spanning trees by Matrix--Tree.
The positive universal assertion is NOT proved by these finite checks.
"""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import product
from math import factorial, prod

LIMIT = 200


def tree_count(vertices, edges):
    """Exact symmetric Gaussian elimination of a Laplacian cofactor.

    Keep only its upper triangle. Products of Schur-complement pivots
    give the determinant; connected input makes every pivot positive.
    The 0-by-0 cofactor has determinant one, including the one-vertex case.
    """
    index = {v: i for i, v in enumerate(vertices)}
    size = len(vertices) - 1
    rows = [defaultdict(Q) for _ in range(size)]
    for u, v in set(tuple(sorted((index[u], index[v]))) for u, v in edges):
        assert u != v
        if u < size:
            rows[u][u] += 1
        if v < size:
            rows[v][v] += 1
        if v < size:
            rows[u][v] -= 1
    det = Q(1)
    for k in range(size):
        pivot = rows[k][k]
        assert pivot > 0, (k, pivot)
        det *= pivot
        tail = [(j, x) for j, x in rows[k].items() if j > k and x]
        tail.sort()
        for pos, (i, x) in enumerate(tail):
            for j, y in tail[pos:]:
                rows[i][j] -= x * y / pivot
        rows[k].clear()
    assert det.denominator == 1
    return det.numerator


def grid(n):
    # Slow coordinates first: for a<=b<=c this keeps elimination bands small.
    vertices = list(product(*(range(s) for s in reversed(n))))
    edges = []
    for v in vertices:
        for i, s in enumerate(reversed(n)):
            if v[i] + 1 < s:
                w = list(v)
                w[i] += 1
                edges.append((v, tuple(w)))
    return vertices, edges


def sigma(n):
    return sum(prod(n[:i] + n[i+1:]) for i in range(len(n)))


def log_majorizes(n, m):
    assert prod(n) == prod(m)
    # Natural log is increasing: compare products of descending prefixes.
    n, m = sorted(n, reverse=True), sorted(m, reverse=True)
    return all(prod(n[:k]) >= prod(m[:k]) for k in range(1, len(n)))


def multiplicity(n):
    return factorial(len(n)) // prod(factorial(n.count(s)) for s in set(n))


def log_interval(x, terms=60):
    """Rational enclosure of log(x), x>=1, by range reduction and atanh.

    For 1<=y<=2, z=(y-1)/(y+1) is at most 1/3 and
    0 <= log(y)-2 sum_{j<M} z^(2j+1)/(2j+1)
      <= 2*z^(2M+1)/((2M+1)*(1-z*z)).
    """
    assert x >= 1
    power = max(0, x.numerator.bit_length() - x.denominator.bit_length())
    if x < 2**power:
        power -= 1
    y = x / 2**power
    assert 1 <= y <= 2

    def unit(y):
        z = (y - 1) / (y + 1)
        low = 2 * sum((z**(2*j+1) / (2*j+1) for j in range(terms)), Q(0))
        rem = 2 * z**(2*terms+1) / ((2*terms+1) * (1-z*z))
        return low, low + rem

    lo, hi = unit(y)
    lo2, hi2 = unit(Q(2))
    return lo + power*lo2, hi + power*hi2


def decimal_enclosure(lo, hi, places=15):
    scale = 10**places
    a, b = (lo*scale).__floor__(), (hi*scale).__ceil__()
    def show(v):
        return f'{v//scale}.{v%scale:0{places}d}'
    return '[' + show(a) + ', ' + show(b) + ']'


def fidelity():
    # Source Section 6: nonstandard periodic quotient vs standard 3x3 torus.
    vs = list(range(9))
    quotient_edges = [(i, (i+s) % 9) for i in vs for s in (1, 3)]
    quotient = tree_count(vs, quotient_edges)
    vs2 = list(product(range(3), repeat=2))
    torus_edges = [(v, tuple((v[j]+(j == i)) % 3 for j in range(2)))
                   for v in vs2 for i in range(2)]
    torus = tree_count(vs2, torus_edges)
    assert (quotient, torus) == (12321, 11664)
    print(f'FIDELITY Section 6 periodic quotient tau={quotient}; C3 box C3 tau={torus}')
    # Independent determinant algorithm on small grids and a nontrivial 24-box.
    import sympy as sp
    for n in [(1,1,1), (1,2,3), (2,2,2), (1,3,4), (2,3,4)]:
        vs, es = grid(n)
        ix = {v: i for i,v in enumerate(vs)}
        lap = sp.zeros(len(vs))
        for u,v in es:
            i,j = ix[u],ix[v]
            lap[i,i] += 1
            lap[j,j] += 1
            lap[i,j] -= 1
            lap[j,i] -= 1
        exact = int(lap[:-1,:-1].det(method='bareiss'))
        assert tree_count(vs, es) == exact
    print('CROSS_CHECK sympy Bareiss cofactors agree on five boxes (N=1,6,8,12,24)')


def main():
    fidelity()
    by_volume = defaultdict(list)
    for a in range(1, LIMIT+1):
        for b in range(a, LIMIT//a+1):
            if a*b*b > LIMIT:
                break
            for c in range(b, LIMIT//(a*b)+1):
                by_volume[a*b*c].append((a,b,c))
    tau = {}
    for n in (n for ns in by_volume.values() for n in ns):
        tau[n] = tree_count(*grid(n))
    assert tau[(1,2,12)] == 2107560
    assert tau[(2,3,4)] == 76752081000
    assert tau[(1,1,64)] == 1
    assert tau[(4,4,4)] == 172685928902844729688524604506636288
    print('CHECK tau(1,2,12)=2107560; tau(2,3,4)=76752081000; sigma deficit=12')
    print('CHECK tau(1,1,64)=1; tau(4,4,4)=172685928902844729688524604506636288; sigma deficit=81')
    pairs = equal = ordered_pairs = 0
    best = None
    for volume, ns in sorted(by_volume.items()):
        for n in ns:
            for m in ns:
                if not log_majorizes(n,m):
                    continue
                deficit = sigma(n)-sigma(m)
                assert deficit >= 0
                if deficit == 0:
                    assert n == m and tau[n] == tau[m]
                    equal += 1
                    continue
                assert tau[m] > tau[n]
                # d=3: c3=log(5/4)/2. Exponentiating gives an INTEGER test.
                assert tau[m]**2 * 4**deficit >= tau[n]**2 * 5**deficit, (n,m)
                pairs += 1
                ordered_pairs += multiplicity(n)*multiplicity(m)
                # Compare log(A)/D to log(B)/E by A**E < B**D.
                if best is None or (tau[m]**best[2] * tau[best[0]]**deficit
                                    < tau[n]**best[2] * tau[best[1]]**deficit):
                    best = n,m,deficit
    print(f'DOMAIN d=3; 1<=N<={LIMIT}; all positive integer triples modulo coordinate permutations')
    print(f'BOXES sorted={len(tau)}; ordered={sum(multiplicity(n) for n in tau)}')
    print(f'PAIRS strict_sorted={pairs}; strict_ordered={ordered_pairs}; equality_sorted={equal}')
    print('INEQUALITY every comparable pair passes tau_m^2*4^deficit >= tau_n^2*5^deficit')
    n,m,deficit = best
    print(f'MIN_RATIO n={n}; m={m}; N={prod(n)}; sigma_deficit={deficit}')
    print(f'MIN_RATIO exact=log({tau[m]}/{tau[n]})/{deficit}')
    lo,hi = log_interval(Q(tau[m],tau[n]))
    lo,hi = lo/deficit,hi/deficit
    print('MIN_RATIO certified_decimal_enclosure='+decimal_enclosure(lo,hi))
    clo,chi = log_interval(Q(5,4))
    print('C3 certified_decimal_enclosure='+decimal_enclosure(clo/2,chi/2))
    assert lo > chi/2
    assert sigma((1,1,1)) == 3 and tau[(1,1,1)] == 1
    print('DEGENERATE all-one equality included; all zero-deficit pairs have equal sorted tuples')
    print('SCOPE finite evidence only; neither a universal proof nor a Lean verification')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
