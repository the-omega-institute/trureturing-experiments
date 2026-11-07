#!/usr/bin/env python3
"""Exact checks for Skoric's Conjecture 1, arXiv:2012.15493v2, equation (35).

The finite checks are evidence about the listed ranges, not a universal proof.
Only integers and fractions.Fraction participate in mathematical comparisons.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import comb


def normalized_means(a):
    s, t = Q(0), Q(0)
    out = [None]
    for k, ak in enumerate(a):
        assert ak > 0
        s += ak
        t += k * ak
        if k:
            out.append(t / (k * s))
    return out


def check_blocks(a):
    """Check the proposed identity, all inward comparisons and block domination."""
    m = len(a) - 1
    for k in range(1, m):
        assert k * a[k] ** 2 >= (k + 1) * a[k - 1] * a[k + 1]
    for r in range(1, m):
        sr, tr = sum(a[:r + 1]), sum(k * a[k] for k in range(r + 1))
        sp, tp = sr + a[r + 1], tr + (r + 1) * a[r + 1]
        d = (r + 1) * sp * tr - r * sr * tp
        negative = (r + 1) * a[r + 1] * sum((r - j) * a[j] for j in range(r))
        assert d == sr * tr - negative
        regrouped = sum(Q(s, 2) * sum(a[i] * a[s - i]
                        for i in range(max(0, s - r), min(r, s) + 1))
                        for s in range(2 * r + 1))
        assert regrouped == sr * tr
        for j in range(r):
            s = r + 1 + j
            endpoint = a[j + 1] * a[r]
            bound = Q(r + 1, j + 1) * a[j] * a[r + 1]
            assert endpoint >= bound
            for i in range(j + 1, r + 1):
                assert a[i] * a[s - i] >= endpoint
            block = Q(s, 2) * sum(a[i] * a[s - i] for i in range(j + 1, r + 1))
            assert block >= (r + 1) * (r - j) * a[j] * a[r + 1]
        assert d >= a[0] * a[1] > 0
    fs = normalized_means(a)
    assert all(fs[r + 1] < fs[r] for r in range(1, m))


# Fidelity: the exact equality in equation (45)
# in the proof of Proposition 1; it does not assume Conjecture 1.
for z in range(1, 101):
    fs = normalized_means([Q(comb(2 * z, k)) for k in range(z + 1)])
    assert fs[z] == 1 / (1 + Q(comb(2 * z, z), 2 ** (2 * z)))
print('FIDELITY: f(2z,z) = 1/(1+binom(2z,z)/2^(2z)), z=1..100; 100 identities')

# Every source-domain pair for all small integer rows, plus the stronger
# all-cap adjacent assertion, including r=1 and the upper endpoint r=n.
pairs, adjacent = 0, 0
for n in range(2, 201):
    a = [Q(comb(n, k)) for k in range(n + 1)]
    fs = normalized_means(a)
    for r in range(1, n):
        assert fs[r + 1] < fs[r]
        adjacent += 1
    for s in range(2, n // 2 + 1):
        for r in range(1, s):
            assert fs[s] < fs[r]
            pairs += 1
    if n <= 12:
        check_blocks(a)
assert normalized_means([Q(comb(4, k)) for k in range(3)])[1:] == [Q(4, 5), Q(8, 11)]
print(f'BINOMIAL: n=2..200; {pairs} source-domain pairs; {adjacent} stronger adjacent checks; all strict')

# Enumerate finite ULC rows via nonincreasing q_k=(k+1)a_(k+1)/a_k.
levels = [Q(1, 3), Q(1, 2), Q(1), Q(3), Q(8)]
rows = 0
for m in range(2, 8):
    for increasing in combinations_with_replacement(levels, m):
        qs = tuple(reversed(increasing))
        a = [Q(1)]
        for k, q in enumerate(qs):
            a.append(a[-1] * q / (k + 1))
        check_blocks(a)
        rows += 1
print(f'ULC: {rows} exact rows, lengths 3..8; identity, inward pairs, block domination and strictness')

# Real-x extension is checked separately; it is not part of the binding claim.
xs = sorted({Q(m, d) for d in (2, 3, 5) for m in range(2 * d, 40 * d + 1)})
real_pairs = 0
for x in xs:
    cap = x.numerator // (2 * x.denominator)
    a = [Q(1)]
    for k in range(cap):
        a.append(a[-1] * (x - k) / (k + 1))
    fs = normalized_means(a)
    for s in range(2, cap + 1):
        for r in range(1, s):
            assert fs[s] < fs[r]
            real_pairs += 1
print(f'REAL-X EXTENSION: {len(xs)} rational x in [2,40]; {real_pairs} source-cap pairs; all strict')

# A separate family checks the proposed structural generalization.
pb_rows = 0
for n in range(2, 8):
    for ps in combinations_with_replacement([Q(1, 5), Q(1, 2), Q(4, 5)], n):
        a = [Q(1)]
        for p in ps:
            b = [Q(0)] * (len(a) + 1)
            for k, ak in enumerate(a):
                b[k] += (1 - p) * ak
                b[k + 1] += p * ak
            a = b
        assert sum(a) == 1
        check_blocks(a)
        pb_rows += 1
print(f'POISSON-BINOMIAL: {pb_rows} exact interior laws, 2..7 Bernoulli coordinates; all strict')
print('SCOPE: finite exact checks passed; no universal claim or Lean proof certified by this program')
print('ALL_CHECKS_OK')
