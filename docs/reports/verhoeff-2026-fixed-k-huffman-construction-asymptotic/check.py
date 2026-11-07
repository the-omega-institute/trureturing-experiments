#!/usr/bin/env python3
"""Exact checks of Verhoeff's balanced, unreduced construction counts.

Source: arXiv:2605.27617v1, Sections 2, 3.1, 3.4 and 8.
Finite checks support the proposed asymptotic; they do not prove its limit.
Only Python's standard library and exact integer/rational arithmetic are used.
"""

from fractions import Fraction as Q
from functools import cache
from heapq import heapify, heappop, heappush
from itertools import combinations_with_replacement
from math import comb, factorial
import json


def huffman(weights):
    """A commutator merges weights a,b to 2(a+b), not a+b."""
    heap = list(weights)
    assert heap and all(w > 0 for w in heap)
    heapify(heap)
    while len(heap) > 1:
        a, b = heappop(heap), heappop(heap)
        heappush(heap, 2 * (a + b))
    return heap[0]


@cache
def optimum(weights):
    """Independent exhaustive optimization over root bipartitions."""
    if len(weights) == 1:
        return weights[0]
    values = []
    # Put the first leaf on the left to remove root-reflection duplicates.
    for mask in range(1, (1 << len(weights)) - 1, 2):
        left = tuple(w for j, w in enumerate(weights) if mask & (1 << j))
        right = tuple(w for j, w in enumerate(weights) if not mask & (1 << j))
        values.append(2 * (optimum(left) + optimum(right)))
    return min(values)


def make_counts(max_k, max_i):
    counts = [[0] * (max_i + 1) for _ in range(max_k + 1)]
    for i in range(max_i + 1):
        n = 2**i
        for k in range(1, min(max_k, n) + 1):
            if k == n:
                value = n
            elif k + 1 == n:
                value = 2 * n
            elif k == 1:
                value = 4**i
            else:
                m = n // 2
                weights = [counts[j][i - 1] + counts[k - j][i - 1]
                           for j in range(max(0, k - m), min(k, m) + 1)]
                value = huffman(weights)
            counts[k][i] = value
    return counts


@cache
def count_by_root_optimization(k, i):
    n = 2**i
    if k == 0 or k > n:
        return 0
    if k == n:
        return n
    if k + 1 == n:
        return 2 * n
    if k == 1:
        return 4**i
    m = n // 2
    weights = tuple(sorted(
        count_by_root_optimization(j, i - 1)
        + count_by_root_optimization(k - j, i - 1)
        for j in range(max(0, k - m), min(k, m) + 1)))
    return optimum(weights)


def coefficients(k):
    """Proposed eventual increment: a[k,i]-6a[k,i-1]=sum alpha[j]a[j,i-1]."""
    alpha = [0] * k
    s = k // 2
    for r in range(1, s):
        alpha[r] += 24 * 4**(r - 1)
        alpha[k - r] += 24 * 4**(r - 1)
    if k % 2 == 0:
        alpha[s] += 2 * 4**s
    else:
        alpha[s] += 4**(s + 1)
        alpha[s + 1] += 4**(s + 1)
    return alpha


def evaluate(poly, x):
    return sum(c * x**j for j, c in enumerate(poly))


def decimal_floor(q, places=6):
    """Display floor(10^places q)/10^places using integers only."""
    assert q >= 0
    scale = 10**places
    v = q.numerator * scale // q.denominator
    return f"{v // scale}.{v % scale:0{places}d}"


def main():
    trials = 0
    for length in range(1, 9):
        for weights in combinations_with_replacement((1, 2, 3, 4, 8), length):
            assert huffman(weights) == optimum(weights), weights
            trials += 1
    assert trials == 1286
    print(f"GREEDY_VS_ALL_ROOT_BIPARTITIONS: {trials} multisets, 1..8 leaves")
    optimum.cache_clear()

    counts = make_counts(12, 160)
    assert counts[2][2] == 32
    assert counts[2][10] == 157048832
    for i in range(2, 161):
        assert counts[2][i] == Q(8, 3) * 6**i - 4 * 4**i
    print("FIDELITY: a[2,2]=32; a[2,10]=157048832; source k=2 formula at 2<=i<=160")

    comparisons = 0
    for k in range(1, 6):
        for i in range(41):
            assert counts[k][i] == count_by_root_optimization(k, i), (k, i)
            comparisons += 1
    print(f"INDEPENDENT_RECURRENCE: {comparisons} pairs, k=1..5, i=0..40")
    for i in range(6, 161):
        assert counts[3][i] == 6**i * (Q(64, 9) * i - Q(7136, 243)) + 24 * 4**i
    print("K3_IDENTITY: 6^i*((64/9)*i-7136/243)+24*4^i, 6<=i<=160")

    p = {1: [Q(0)], 2: [Q(8, 3)]}
    b = {1: Q(1), 2: Q(-4)}
    polynomial_results = []
    for k in range(3, 13):
        alpha = coefficients(k)
        for i in range(1, 161):
            if 2**(i - 1) < k:
                continue
            total_lower = 2 * sum(counts[j][i - 1] for j in range(1, k))
            increment = counts[k][i] - 6 * counts[k][i - 1]
            assert 2 * total_lower <= increment <= 2**k * total_lower
        matches = [i for i in range(1, 161)
                   if counts[k][i] - 6 * counts[k][i - 1]
                   == sum(alpha[j] * counts[j][i - 1] for j in range(1, k))]
        assert all(i in matches for i in range(40, 161)), k
        start = 160
        while start - 1 in matches:
            start -= 1
        b[k] = -sum(alpha[j] * b[j] for j in range(1, k)) / 2
        degree = k - 2
        rhs = [Q(0)] * degree
        for j in range(1, k):
            for t, c in enumerate(p[j]):
                for h in range(t + 1):
                    rhs[h] += Q(alpha[j], 6) * c * comb(t, h) * (-1)**(t - h)
        poly = [Q(0)] * (degree + 1)
        for h in range(degree - 1, -1, -1):
            known = sum(poly[t] * comb(t, h) * (-1)**(t - h + 1)
                        for t in range(h + 2, degree + 1))
            poly[h + 1] = (rhs[h] - known) / (h + 1)
        # Match the free integration constant at one index; check other indices separately.
        poly[0] = Q(counts[k][40], 6**40) - b[k] * Q(2, 3)**40 - evaluate(poly, 40)
        expected = Q(4**k, 9 * factorial(k - 2))
        assert poly[-1] == expected
        for i in range(40, 161):
            assert counts[k][i] == 6**i * evaluate(poly, i) + b[k] * 4**i
        p[k] = poly
        polynomial_results.append({"k": k, "leading_coefficient": str(poly[-1]),
                                   "b": str(b[k]), "P_low_to_high": list(map(str, poly)),
                                   "increment_matches_through_160_from": start})
    print("EVENTUAL_FORM_FINITE_CHECK: k=3..12, 40<=i<=160; exact leading coefficients")
    print("POLYNOMIALS " + json.dumps(polynomial_results, sort_keys=True))

    ratios = []
    for k in (3, 4, 5):
        c = Q(4**k, 9 * factorial(k - 2))
        for i in range(41):
            if k > 2**i:
                continue
            ratio = Q(counts[k][i], 1) / (c * 6**i * i**(k - 2))
            ratios.append({"k": k, "i": i, "a": counts[k][i],
                           "ratio_exact": str(ratio), "ratio_decimal_floor": decimal_floor(ratio)})
        sampled = [Q(counts[k][i], 1) / (c * 6**i * i**(k - 2)) for i in (10, 20, 40, 80, 160)]
        assert all(0 < x < 1 for x in sampled)
        assert all(x < y for x, y in zip(sampled, sampled[1:]))
        print(f"RATIO_TAIL k={k}: " + "; ".join(
            f"i={i}: {q} (floor6={decimal_floor(q)})"
            for i, q in zip((10, 20, 40, 80, 160), sampled)))
    print("COUNTS_AND_RATIOS_I_LE_40 " + json.dumps(ratios, sort_keys=True))
    print("SCOPE: finite exact checks only; no proof of the universal limit or stabilization")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
