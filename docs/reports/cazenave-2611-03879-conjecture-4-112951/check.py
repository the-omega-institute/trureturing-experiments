#!/usr/bin/env python3
"""Exact finite checks of Cazenave arXiv:2611.03879v1, Conjecture 4.

All certificates use integers and fractions. SymPy computes the source's K2
spectra exactly. No numerical eigenvalues or floating-point tolerances are used.
The general mathematical argument is described in README.md; this script checks
its three sufficient certificates separately on each enumerated graph.
"""
from collections import Counter
from fractions import Fraction
from functools import cache
import gzip
import hashlib
from math import isqrt
from pathlib import Path
import sympy as sp

EXPECTED = {2: 1, 3: 2, 4: 6, 5: 21, 6: 112, 7: 853, 8: 11117, 9: 261080}
HERE = Path(__file__).resolve().parent


@cache
def square_part(k):
    """Return t,s with k=t*t*s and s squarefree."""
    assert k > 0
    t, s, p = 1, k, 2
    while p * p <= s:
        while s % (p * p) == 0:
            t *= p
            s //= p * p
        p += 1
    assert t * t * s == k
    return t, s


@cache
def radical_interval(s, bits):
    """Certified rational interval containing sqrt(s)."""
    q = 1 << bits
    a = isqrt(s * q * q)
    lo = Fraction(a, q)
    hi = lo if a * a == s * q * q else Fraction(a + 1, q)
    assert lo * lo <= s <= hi * hi
    return lo, hi


@cache
def randic_certificate(products, delta):
    """Certify sum_edges 1/sqrt(d_i*d_j) - sqrt(delta) >= 0."""
    coeff = Counter()
    for k, multiplicity in products:
        t, s = square_part(k)
        coeff[s] += Fraction(multiplicity, t * s)
    t, s = square_part(delta)
    coeff[s] -= t
    coeff = {s: c for s, c in coeff.items() if c}
    if not coeff:
        return 'equal'
    for bits in (8, 16, 32, 64, 128, 256):
        lower = Fraction(0)
        for s, c in coeff.items():
            lo, hi = radical_interval(s, bits)
            lower += c * (lo if c > 0 else hi)
        if lower > 0:
            return 'positive'
    raise AssertionError(('uncertified Randic bound', products, delta, coeff))


def qnonneg(a, b, delta):
    """Decide a+b*sqrt(delta)>=0 by integer sign and square tests."""
    if b >= 0:
        return a >= 0 or b * b * delta >= a * a
    return a >= 0 and a * a >= b * b * delta


def graph6(raw):
    vals = [c - 63 for c in raw]
    n = vals[0]
    assert 2 <= n <= 9
    bits = ((v >> k) & 1 for v in vals[1:] for k in range(5, -1, -1))
    edges = [(i, j) for j in range(1, n) for i in range(j) if next(bits)]
    neighbors = [set() for _ in range(n)]
    for i, j in edges:
        neighbors[i].add(j)
        neighbors[j].add(i)
    reached, pending = {0}, [0]
    while pending:
        i = pending.pop()
        for j in neighbors[i] - reached:
            reached.add(j)
            pending.append(j)
    assert len(reached) == n, 'disconnected input'
    return n, edges, neighbors


def certificates(edges, neighbors):
    """Return three exact certificates implying the target via spectral facts.

    1. R>=sqrt(delta), with grouped radicals and certified rational intervals.
    2. A*h <= c*h for h_u=sqrt(delta), h_i=1 otherwise,
       c=sqrt(delta)+d-1. A nonnegative matrix has spectral radius <=c:
       similarity by diag(h) gives a nonnegative matrix with row sums <=c.
    3. The u,v principal block of L-(d-1)I is PSD. Interlacing gives
       the second largest eigenvalue of L >=d-1.
    Together these imply lambda_1 <= R+mu_2 on THIS graph.
    """
    n = len(neighbors)
    degrees = list(map(len, neighbors))
    assert min(degrees) > 0
    u = max(range(n), key=degrees.__getitem__)
    v = max((i for i in range(n) if i != u), key=degrees.__getitem__)
    delta, d = degrees[u], degrees[v]
    assert delta >= d >= 1
    products = tuple(sorted(Counter(degrees[i] * degrees[j] for i, j in edges).items()))
    rcert = randic_certificate(products, delta)
    # Residual c*h-A*h is represented as a+b*sqrt(delta).
    for i in range(n):
        if i == u:
            # c*h_u = delta+(d-1)*sqrt(delta); (A*h)_u=delta.
            a, b = 0, d - 1
        else:
            adjacent = int(u in neighbors[i])
            # (A*h)_i=degree_i-adjacent+adjacent*sqrt(delta).
            a, b = d - 1 - degrees[i] + adjacent, 1 - adjacent
        assert qnonneg(a, b, delta), ('adjacency certificate', i, a, b, delta)
    # Principal block [[delta-d+1,-epsilon],[-epsilon,1]].
    epsilon = int(v in neighbors[u])
    top, bottom, off = delta - d + 1, 1, -epsilon
    assert top >= 0 and bottom >= 0 and top * bottom - off * off >= 0
    return rcert


def fidelity():
    # Section 7: "Equality holds for K_2."
    adjacency = sp.Matrix([[0, 1], [1, 0]])
    laplacian = sp.diag(1, 1) - adjacency
    ae = sorted(adjacency.eigenvals().keys(), reverse=True)
    le = sorted(laplacian.eigenvals().keys(), reverse=True)
    randic = 1 / sp.sqrt(sp.Integer(1) * sp.Integer(1))
    slack = sp.simplify(randic + le[1] - ae[0])
    assert ae == [1, -1] and le == [2, 0] and randic == 1 and slack == 0
    print('FIDELITY K2: adjacency=[1,-1] Laplacian=[2,0] R=1 slack=0')
    # Section 1, classical fact F6: R(K_n)=n/2.
    for n in range(2, 10):
        r = Fraction(n * (n - 1), 2) / (n - 1)
        assert r == Fraction(n, 2)
    print('FIDELITY K_n: R=n/2 for 2<=n<=9 (exact fractions)')


def main():
    fidelity()
    counts, hashes, strict, equal = Counter(), {}, Counter(), Counter()
    streamhash = hashlib.sha256()
    with gzip.open(HERE / 'connected-2-9.g6.gz', 'rb') as stream:
        for line in stream:
            raw = line.rstrip(b'\n')
            n, edges, neighbors = graph6(raw)
            counts[n] += 1
            hashes.setdefault(n, hashlib.sha256()).update(line)
            streamhash.update(line)
            rc = certificates(edges, neighbors)
            (strict if rc == 'positive' else equal)[n] += 1
    assert dict(counts) == EXPECTED, (counts, EXPECTED)
    for n in EXPECTED:
        print(f'ORDER {n}: connected_unlabeled={counts[n]} target_certificates=PASS '
              f'R_gt_sqrt_delta={strict[n]} R_eq_sqrt_delta={equal[n]} '
              f'graph6_sha256={hashes[n].hexdigest()}')
    print(f'TOTAL_CONNECTED_UNLABELED={sum(counts.values())}')
    print(f'GRAPH6_STREAM_SHA256={streamhash.hexdigest()}')
    print('MINIMUM_TARGET_SLACK=0 (all graphs certified >=0; exact equality at K2)')
    print('EVIDENCE=exact certificates plus classical spectral comparison/interlacing; no Lean proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
