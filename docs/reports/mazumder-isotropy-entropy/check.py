#!/usr/bin/env python3
"""Exact checks for Mazumder, arXiv:2610.01024v2, Open Problem 9.

No numerical logarithms are used. Entropy ceilings are decided by rational
powers: ceil(-log_2(m)/log_2(4/3))=e iff (3/4)^e <= m < (3/4)^(e-1).
The endpoint inequalities used below are strict.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations
from time import perf_counter_ns

MASK = int('145f6c4172ac2b1623b67fb656ecd3e07d5cd1e7a6fca87bd0498a9b82355a721ffeb93e3539cc64b087c3e09732045fe1d01b842becc838e1949d2b3c5c0e7f01d304efb507a16262e07f9db8964526e67edcc0dbc2009ac06189cfd6a9f09b49dff3e9d3d52fc2358fc98ddc7bfa79f2a8aaab616cb98521e', 16)
EXPECTED = {0: 1, 12: 31, 14: 10970, 16: 292882, 18: 220404}


def triples(n, mask):
    return [t for r, t in enumerate(combinations(range(n), 3)) if (mask >> r) & 1]


def elementary_slices(n, ts):
    slices = [[0] * n for _ in range(n)]
    for i, j, k in ts:
        for a, b, c in [(i, j, k), (j, i, k), (k, i, j)]:
            slices[a][b] ^= 1 << c
            slices[a][c] ^= 1 << b
    return slices


def rank_insert(rows):
    basis = [0] * len(rows)
    rank = 0
    for v in rows:
        while v:
            p = v.bit_length() - 1
            if basis[p]:
                v ^= basis[p]
            else:
                basis[p] = v
                rank += 1
                break
    return rank


def histogram_gray(n, slices):
    h = Counter({0: 1})
    rows = [0] * n
    for step in range(1, 1 << n):
        i = (step & -step).bit_length() - 1
        rows = [x ^ y for x, y in zip(rows, slices[i])]
        h[rank_insert(rows)] += 1
    return dict(sorted(h.items()))


def histogram_independent(n, mask):
    # Reconstruct directly with three nested loops, rather than combinations.
    # For each output row build the list of coefficient-bearing direction bits.
    coefficients = [[0] * n for _ in range(n)]
    r = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if (mask >> r) & 1:
                    coefficients[j][k] ^= 1 << i
                    coefficients[k][j] ^= 1 << i
                    coefficients[i][k] ^= 1 << j
                    coefficients[k][i] ^= 1 << j
                    coefficients[i][j] ^= 1 << k
                    coefficients[j][i] ^= 1 << k
                r += 1
    # Lookup tables are exact parity evaluations, with a different traversal.
    lo = min(n, 10)
    tables = []
    for row in coefficients:
        tables.append((
            [sum((((a & c).bit_count() & 1) << j) for j, c in enumerate(row))
             for a in range(1 << lo)],
            [sum(((((a << lo) & c).bit_count() & 1) << j) for j, c in enumerate(row))
             for a in range(1 << (n - lo))]))
    h = Counter()
    for a in range(1 << n):
        rows = [l[a & ((1 << lo) - 1)] ^ u[a >> lo] for l, u in tables]
        rank = 0
        # Lowest-column row pivoting, independent of highest-bit insertion.
        for c in range(n):
            p = next((p for p in range(rank, n) if (rows[p] >> c) & 1), None)
            if p is None:
                continue
            rows[rank], rows[p] = rows[p], rows[rank]
            v = rows[rank]
            for p in range(rank + 1, n):
                if (rows[p] >> c) & 1:
                    rows[p] ^= v
            rank += 1
        h[rank] += 1
    return dict(sorted(h.items()))


def slice_sum(n, ts):
    return sum((Fraction(count, 1 << rank)
                for rank, count in histogram_gray(n, elementary_slices(n, ts)).items()), Fraction())


def phase_values(n, ts):
    supports = [sum(1 << i for i in t) for t in ts]
    return [-1 if sum((x & t) == t for t in supports) % 2 else 1
            for x in range(1 << n)]


def pauli_moment_phase(n, ts):
    signs = phase_values(n, ts)
    total = 0
    for a in range(1 << n):
        for b in range(1 << n):
            # X^a Z^b differs from a Hermitian Pauli only by a unit phase;
            # its fourth absolute moment is unchanged.
            w = sum(signs[x] * signs[x ^ a] * (-1 if (b & x).bit_count() % 2 else 1)
                    for x in range(1 << n))
            total += w ** 4
    return Fraction(total, 1 << (5 * n))


def pauli_moment_choi(n, ts):
    # Direct 2n-qubit Pauli sum, using the sparse support |x>|x>.
    signs = phase_values(n, ts)
    size = 1 << n
    total = 0
    support = {x | (x << n): signs[x] for x in range(size)}
    for a in range(1 << (2 * n)):
        for b in range(1 << (2 * n)):
            w = sum(s * support.get(x ^ a, 0) * (-1 if (b & x).bit_count() % 2 else 1)
                    for x, s in support.items())
            total += w ** 4
    return Fraction(total, (1 << (2 * n)) * size ** 4)


def strict_ceiling(m, e):
    assert Fraction(3, 4) ** e < m < Fraction(3, 4) ** (e - 1)


def main():
    # Fidelity: Appendix F states M2(CCZ)=5-log2(11), E=4, F=7.
    ccz = [(0, 1, 2)]
    ccz_m = pauli_moment_phase(3, ccz)
    ccz_j = pauli_moment_choi(3, ccz)
    assert ccz_m == ccz_j == Fraction(11, 32)
    assert histogram_gray(3, elementary_slices(3, ccz)) == {0: 1, 2: 7}
    strict_ceiling(ccz_j, 4)
    print('FIDELITY: direct Choi and phase Pauli sums m2(CCZ)=11/32; M2=5-log2(11); E=4; F=7', flush=True)
    # Exhaust every homogeneous cubic through n=4, including zero.
    checked = 0
    for n in range(5):
        count = len(list(combinations(range(n), 3)))
        for mask in range(1 << count):
            ts = triples(n, mask)
            assert pauli_moment_phase(n, ts) == slice_sum(n, ts) / (1 << n)
            checked += 1
    print(f'FOURTH-MOMENT IDENTITY: all {checked} cubics on n=0,...,4 agree exactly', flush=True)
    # Direct-sum formula on every direction of a 7-variable small example.
    t4 = triples(4, 11)
    t7 = t4 + [tuple(i + 4 for i in t) for t in ccz]
    assert slice_sum(7, t7) == slice_sum(4, t4) * slice_sum(3, ccz)
    print('DIRECT SUM: complete 7-variable slice sum factors exactly', flush=True)
    ts = triples(19, MASK)
    assert len(ts) == 499 and MASK.bit_length() <= 969
    start = perf_counter_ns()
    h1 = histogram_gray(19, elementary_slices(19, ts))
    print(f'WITNESS: Gray/highest-pivot histogram={h1}; elapsed_ns={perf_counter_ns()-start}', flush=True)
    start = perf_counter_ns()
    h2 = histogram_independent(19, MASK)
    print(f'WITNESS: binary/lowest-column histogram={h2}; elapsed_ns={perf_counter_ns()-start}', flush=True)
    assert h1 == h2 == EXPECTED
    assert sum(h1.values()) == 1 << 19
    s = sum((Fraction(v, 1 << r) for r, v in h1.items()), Fraction())
    assert s == Fraction(457895, 65536) < 7
    m = s / (1 << 19)
    assert m == Fraction(457895, 1 << 35)
    strict_ceiling(m, 40)
    assert 3 ** 39 == 4052555153018976267 > 7 * 2 ** 59 == 4035225266123964416
    assert 3 ** 40 == 12157665459056928801 < 457895 * 2 ** 45 == 16110748057615728640
    print(f'EXACT: radical cardinality=1; d=0; F=39; S={s}; Choi m2={m}; E=40', flush=True)
    # General family argument uses m_k=m^k, F_k=38k+1.
    # m<(3/4)^39 implies m^k<(3/4)^(39k)<=(3/4)^(38k+1).
    assert m < Fraction(3, 4) ** 39
    for k in range(1, 33):
        assert 39 * k >= 38 * k + 1
        assert m ** k < Fraction(3, 4) ** (38 * k + 1)
    print('FAMILY: exact base inequality m<(3/4)^39; direct sums multiply moments; 39k>=38k+1 for every k>=1 (power checks k=1,...,32)', flush=True)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
