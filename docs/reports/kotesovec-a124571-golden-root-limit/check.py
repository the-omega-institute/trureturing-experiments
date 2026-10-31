#!/usr/bin/env python3
"""Exact finite checks for Hanna--Kotesovec, OEIS A124571 revision 20.

These checks test source fidelity and the proposed all-index proof's identities.
They do not establish the infinite limit. No floating point is used.
"""
from fractions import Fraction
from math import comb

PREFIX = [
    1, 1, 1, 3, 4, 6, 12, 21, 33, 53, 90, 151, 253, 426, 701, 1151,
    1900, 3123, 5162, 8553, 14092, 23223, 38296, 62963, 103458, 170056,
    279140, 457833, 751033, 1231671, 2019090, 3309710, 5424315,
    8886249, 14553015, 23826952, 38997232, 63806394,
]
LIMIT = 512


def mul(a, b, cap=None):
    degree = len(a) + len(b) - 2
    if cap is not None:
        degree = min(degree, cap)
    out = [0] * (degree + 1)
    for i, ai in enumerate(a[:degree + 1]):
        if ai:
            for j in range(min(len(b) - 1, degree - i) + 1):
                out[i + j] += ai * b[j]
    return out


def add(a, b):
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
            for i in range(max(len(a), len(b)))]


def coeff(a, i):
    return a[i] if 0 <= i < len(a) else 0


def power(a, m):
    out = [1]
    for _ in range(m):
        out = mul(out, a)
    return out


def qsqrt5_eval_at_rho(p):
    """Return integers (a,b) for p(rho)=a+b*rho, rho^2=1-rho."""
    a, b = 0, 0
    for c in reversed(p):
        a, b = b + c, a - b
    return a, b


def tree_shape(m):
    """An explicit balanced unary/binary tree: (leaves, internal nodes, depth)."""
    t = (m - 1).bit_length()

    def build(leaves, depth):
        if depth == 0:
            assert leaves == 1
            return 1, 0, 0
        capacity = 1 << (depth - 1)
        if leaves <= capacity:
            count, internal, height = build(leaves, depth - 1)
            return count, internal + 1, height + 1
        left = build(capacity, depth - 1)
        right = build(leaves - capacity, depth - 1)
        assert left[2] == right[2]
        return left[0] + right[0], left[1] + right[1] + 1, left[2] + 1

    return build(m, t)


def nm(m, k, N):
    return m * (k + 1) - 1 - (m - 1).bit_count() - N


def nu2(m):
    assert m >= 1
    return (m & -m).bit_length() - 1


def root_bracket(a, n, denominator=10**10):
    """Consecutive rational endpoints, checked by cross-multiplied powers."""
    target = a * denominator**n
    lo, hi = 0, 2 * denominator
    while hi**n <= target:
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid**n <= target:
            lo = mid
        else:
            hi = mid
    assert lo**n <= target < hi**n
    return Fraction(lo, denominator), Fraction(hi, denominator)


def main():
    # Source FORMULA, independently evaluated in a truncated integer series ring.
    S = [0, 1] + [0] * (LIMIT - 1)
    snapshots = [S[:]]
    for N in range(LIMIT):
        square = mul(S, S, LIMIT - N)
        new = S[:]
        for j, value in enumerate(square):
            new[N + j] += value
        assert new[:min(N + 2, LIMIT + 1)] == S[:min(N + 2, LIMIT + 1)]
        S = new
        if N < 8:
            snapshots.append(S[:])
    assert S[1:len(PREFIX) + 1] == PREFIX
    assert all(c > 0 for c in S[1:])
    print(f'FIDELITY: all {len(PREFIX)} published terms reproduced exactly')
    print(f'STABILIZATION: source updates preserve degrees < N+2, N=0..{LIMIT-1}')

    # Literal source b(1)=x^2, b(n+1)=b(n)^2+x*b(n), independent of S.
    b = [0, 0, 1]
    P = [[0, 1]]
    for N in range(9):
        p = P[N]
        assert b == [0] + p
        for k in range(min(LIMIT - N, len(p))):
            assert coeff(p, N + k) == snapshots[N][k]
        D = 2**(N + 1) - 1
        assert len(p) - 1 == D and p[-1] == 1
        assert all(coeff(p, k) > 0 for k in range(N + 1, D + 1))
        assert all(coeff(p, k) == 0 for k in range(N + 1))
        assert qsqrt5_eval_at_rho(p) == (0, 1)
        if N < 8:
            P.append([0] + add(p, mul(p, p)))
            b = add(mul(b, b), [0] + b)
    print('NORMALIZATION: b(N+1)=x*P_N and P_N=x^N*S_N, N=0..8')
    print('CRITICAL IDENTITY: P_N(rho)=rho in Z[rho]/(rho^2+rho-1), N=0..8')
    print('SUPPORT: monic, degree 2^(N+1)-1, positive throughout [N+1,D_N], N=0..8')

    tree_cases = 0
    for m in range(1, 4097):
        leaves, v, t = tree_shape(m)
        assert leaves == m and t == (m - 1).bit_length()
        assert v == m - 1 + t - (m - 1).bit_count()
        assert v == sum((m + (1 << j) - 1) // (1 << j) for j in range(1, t + 1))
        for k in range(2, 21):
            assert nm(m + 1, k, 0) - nm(m, k, 0) == k + nu2(m)
        tree_cases += 1
    print(f'TREES/GAPS: {tree_cases} explicit trees; exact digit/gap identities for k=2..20')

    dominance_cases = 0
    lower_cases = 0
    for N in range(1, 4):
        p = P[N]
        for m in range(1, 25):
            _, v, t = tree_shape(m)
            witness = [0] * v + power(p, m)
            target = P[N + t]
            assert len(witness) <= len(target)
            assert all(c <= coeff(target, j) for j, c in enumerate(witness))
            dominance_cases += 1
            for k in range(N + 1, len(p) - 1):
                c, beta = p[k], p[k + 1]
                for d in range(m + 1):
                    index = nm(m, k, N) + d
                    assert 1 <= index <= LIMIT
                    bound = comb(m, d) * c**(m-d) * beta**d
                    assert S[index] >= bound >= c**(m-d)
                    lower_cases += 1
    print(f'TREE DOMINANCE: {dominance_cases} coefficientwise polynomial inequalities')
    print(f'ALL-INDEX LOWER-BOUND INGREDIENT: {lower_cases} adjacent-coefficient inequalities')

    for n in [64, 128, 256, 512]:
        lo, hi = root_bracket(S[n], n)
        assert Fraction(8, 5) <= lo < hi <= Fraction(17, 10)
        print(f'ROOT BRACKET n={n}: {lo} <= a(n)^(1/n) < {hi} (integer powers)')
    print('FINITE CHECKS ONLY: the ordinary infinite limit still requires the analytic proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
