#!/usr/bin/env python3
"""Exact checks for Hanna's OEIS A386648, revision 9.

Only Python integers are used. A list f represents n! [x^n] F(x).
The implicit solution is computed from the original NAME/FORMULA (1),
without using the proposed closed congruence to generate its coefficients.
OEIS source: https://oeis.org/A386648 (CC BY-SA 4.0).
Verification code: Apache-2.0, under the experiments repository's LICENSE.
"""
from math import comb, factorial, isqrt

N = 36
SOURCE_20 = [
    1, 1, 10, 63, 806, 9485, 161540, 2752155, 59021506, 1310350929,
    34127883032, 934203381287, 28851653188814, 942891341070237,
    33962521081521076, 1297690184864525043, 53814377103189792794,
    2366017294084046632937, 111607600081334119137488,
    5565256312162642805357247,
]
CHOOSE = [[comb(n, i) for i in range(n + 1)] for n in range(N + 1)]


def unit(n):
    return [1] + [0] * n


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def mul(a, b, n):
    return [sum(CHOOSE[d][i] * a[i] * b[d - i]
                for i in range(d + 1) if a[i] and b[d - i])
            for d in range(n + 1)]


def power(a, k, n):
    result = unit(n)
    while k:
        if k % 2:
            result = mul(result, a, n)
        k //= 2
        if k:
            a = mul(a, a, n)
    return result


def divide_exact(a, d):
    assert d > 0
    assert all(x % d == 0 for x in a), ("nonintegral divided power", d)
    return [x // d for x in a]


def exp_series(a, n):
    assert a[0] == 0
    result, term = unit(n), unit(n)
    for j in range(1, n + 1):
        term = mul(term, a, n)
        result = add(result, divide_exact(term, factorial(j)))
    return result


def w_series(n):
    # Source EXAMPLE, RELATED SERIES: -LambertW(-x).
    return [0] + [-k ** (k - 1) for k in range(1, n + 1)]


def original_sum(a, n):
    """Sum_m (A^m + W)^m/m!, truncated after degree n."""
    w = w_series(n)
    result, a_to_m = unit(n), unit(n)
    for m in range(1, n + 1):
        a_to_m = mul(a_to_m, a, n)
        term = power(add(a_to_m, w), m, n)
        result = add(result, divide_exact(term, factorial(m)))
    return result


def implicit_solution(n):
    a = [0] * (n + 1)
    for d in range(1, n + 1):
        # The m=1 term has coefficient a[d]+w[d]. All other terms
        # depend on strictly smaller positive-degree coefficients of A.
        a[d] = -original_sum(a[:d + 1], d)[d]
    return a


def main():
    a = implicit_solution(N)
    assert a[1:21] == SOURCE_20
    print("FIDELITY: all 20 displayed OEIS A386648 values match")
    assert original_sum(a, N) == unit(N)
    print(f"ORIGINAL_EQUATION: coefficients 0..{N} equal 1,0,...,0 exactly")
    print(f"INTEGRALITY: all {N} computed EGF coefficients are integers")

    w = w_series(N)
    e = exp_series(w, N)
    x = [0, 1] + [0] * (N - 1)
    minus_x = [-v for v in x]
    assert mul(w, e, N) == minus_x
    print(f"LAMBERT_FIDELITY: W*exp(W)=-x through degree {N}")

    tail = [0] * (N + 1)
    formula2 = [0] * (N + 1)
    for k in range(isqrt(N) + 1):
        front = divide_exact(power(a, k * k, N), factorial(k))
        term = mul(front, exp_series(mul(w, power(a, k, N), N), N), N)
        formula2 = add(formula2, term)
        if k >= 2:
            tail = add(tail, term)
    assert formula2 == unit(N)
    assert all(t % 6 == 0 for t in tail)
    print(f"FORMULA_2: equivalent source equation verified through degree {N}")
    print(f"TAIL_REDUCTION: every EGF tail coefficient is divisible by 6 through degree {N}")

    # Predicted identity: A = 1-W-2x+(x-1)*exp(W) modulo 6 in the EGF lattice.
    one = unit(N)
    x_minus_one = [u - v for u, v in zip(x, one)]
    last = mul(x_minus_one, e, N)
    predicted = [one[d] - w[d] - 2 * x[d] + last[d] for d in range(N + 1)]
    assert all((a[d] - predicted[d]) % 6 == 0 for d in range(N + 1))
    for d in range(2, N + 1):
        closed = (d - 1) ** (d - 1) + d ** (d - 1) - d * (d - 2) ** (d - 2)
        assert (a[d] - closed) % 6 == 0
    print(f"CLOSED_CONGRUENCE: coefficient formula holds for n=2..{N}")
    pattern = [2, 3, 4, 3, 2, 5]
    assert all(a[d] % 6 == pattern[(d - 7) % 6] for d in range(7, N + 1))
    print(f"CONJECTURE: {N - 6} residue checks at n=7..{N}; mismatches=0")
    print("RESIDUES_7_TO_36:", [a[d] % 6 for d in range(7, N + 1)])
    # This only tests the elementary modular formula; it does not extend
    # the independently computed implicit-series range beyond N.
    for d in range(7, 10001):
        closed_mod = (pow(d - 1, d - 1, 6) + pow(d, d - 1, 6)
                      - d * pow(d - 2, d - 2, 6)) % 6
        assert closed_mod == pattern[(d - 7) % 6]
    print("ARITHMETIC_HELPER: closed formula has the predicted residues for n=7..10000")
    print("BOUNDARY: finite exact checks only; no all-index Lean proof is claimed")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
