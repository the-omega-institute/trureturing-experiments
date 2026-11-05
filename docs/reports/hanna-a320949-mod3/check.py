#!/usr/bin/env python3
"""Exact finite checks for Hanna's OEIS A320949; no formal-proof claim.

Source: https://oeis.org/A320949 (revision 12) and its b-file, n=0..24.
Arithmetic uses Python integers and fractions.Fraction exclusively.
"""
from fractions import Fraction as F
from math import comb

N = 24
SOURCE = [1, 3, 21, 360, 9957, 374253, 17525043, 974215107, 62408814861, 4516276087767, 363954839947020, 32309808570484494, 3132592962836443359, 329373550369411338717, 37333941689722714789794, 4538664786712251729910968, 589148225684126151109083696, 81336502478379195711058521051, 11901242590779904151733826269744, 1839862821593669244253225972075464, 299668244116570657639682105216317053, 51292095416570172123936037156127835531, 9204662738338482512852085493152759903051, 1728201764798799473710128066622726184438841, 338819858653227423213915285075931670188366500]


def newton(values):
    out = []
    while values:
        out.append(values[0])
        values = [b - a for a, b in zip(values, values[1:])]
    return out


def poly_mul(a, b):
    c = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return c


def series_mul(a, b, degree):
    return [sum(a[i] * b[n - i] for i in range(n + 1))
            for n in range(degree + 1)]


def residual(b, h, degree):
    powers = [[1] + [0] * degree]
    for _ in range(degree):
        powers.append(series_mul(powers[-1], b, degree))
    return [sum(h[k][r] * powers[r][n - k]
                for k in range(n + 1) for r in range(n - k + 1))
            - (n == 0) for n in range(degree + 1)]


def solve(h, degree):
    b = [0] * (degree + 1)
    for n in range(1, degree + 1):
        b[n] = residual(b[:n + 1], h, n)[n]
    assert residual(b, h, degree) == [0] * (degree + 1)
    return b


def main():
    # Newton route: P_k(t) = binom(t^2,k), deg P_k = 2k.
    # H_k(B) = sum_m c[k,m] (1+B)^(-m-1).
    c = [newton([comb(j * j, k) for j in range(2 * k + 1)])
         for k in range(N + 1)]
    h = [[(-1) ** r * sum(cm * comb(m + r, r)
                          for m, cm in enumerate(c[k]))
          for r in range(N - k + 1)] for k in range(N + 1)]

    # Independent rational-polynomial route to the original geometric sum.
    # M_d = sum_{j>=0} j^d/2^(j+1); shifting the index gives this recursion.
    moments = [1]
    for d in range(1, 2 * N + 1):
        moments.append(sum(comb(d, i) * moments[i] for i in range(d)))
    p = [[F(1)]]
    rising = [[F(1)]]
    for k in range(1, N + 1):
        p.append([a / k for a in poly_mul(p[-1], [F(-(k - 1)), F(0), F(1)])])
        rising.append([a / k for a in poly_mul(rising[-1], [F(k), F(1)])])
    hh = []
    for k in range(N + 1):
        row = []
        assert all(p[k][i] == 0 for i in range(1, len(p[k]), 2))
        for r in range(N - k + 1):
            product = poly_mul(p[k], rising[r])
            exact_sum = sum(a * moments[d] for d, a in enumerate(product))
            value = F((-1) ** r, 2 ** r) * exact_sum
            assert value == h[k][r], (k, r, value, h[k][r])
            row.append(value)
        hh.append(row)
    print(f"Original geometric-moment and Newton transforms agree for all {sum(map(len, h))} pairs k+r <= {N}.")

    b = solve(h, N)
    bb = solve(hh, N)
    assert b == bb
    a = [1] + b[1:]
    assert a == SOURCE
    assert all(isinstance(v, int) and v % 3 == 0 for v in a[1:])
    print(f"Fidelity: a(0)..a({N}) exactly match the OEIS b-file.")
    print(f"Settlement check: all {N} positive-degree coefficients are integers divisible by 3.")
    print("coefficients =", a)

    # A second check of the triangular dependence, on arbitrary integral prefixes.
    trial = [0] + [(-1) ** n * (n + 2) for n in range(1, N + 1)]
    for n in range(1, N + 1):
        old = residual(trial[:n + 1], h, n)[n]
        changed = trial[:n + 1]
        changed[n] += 1
        assert residual(changed, h, n)[n] == old - 1
    print(f"Triangular coefficient slope -1 verified for n=1..{N} on a separate integral prefix.")

    # Exact reciprocity for the integer-valued even P_k, including k! divisible by 3.
    for k in range(61):
        ck = newton([comb(j * j, k) for j in range(2 * k + 1)])
        pzero = int(k == 0)
        total = sum(ck)
        assert total % 3 == pzero
        for z in (F(2), F(3), F(-2), F(2, 3)):
            def q(w):
                return sum(F(cm) * w ** m / (1 - w) ** (m + 1)
                           for m, cm in enumerate(ck))
            assert q(z) + q(1 / z) == pzero, (k, z)
        assert sum(cm * (-1) ** (m + 1) * 2 ** m
                   for m, cm in enumerate(ck)) + 2 * total == pzero
    print("Reciprocity and the mod-3 Newton sum checked for P_k, k=0..60, at four exact rational points.")
    print("No infinite-range conclusion is inferred from these finite checks.")
    print("ALL_CHECKS_OK")


if __name__ == '__main__':
    main()
