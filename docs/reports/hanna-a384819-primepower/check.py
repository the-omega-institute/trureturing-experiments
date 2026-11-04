#!/usr/bin/env python3
"""Exact, offline checks of Hanna's OEIS A384819 prime-power conjecture.

A384819 source: https://oeis.org/A384819, revision 14.
Published table: https://oeis.org/A384819/b384819.txt (all 520 terms).
A384820 source: https://oeis.org/A384820, revision 12 (41 displayed terms).
OEIS data are used under the OEIS End-User License Agreement:
https://oeis.org/LICENSE . Attribution: Paul D. Hanna and The OEIS Foundation Inc.
The computations are finite checks, not an all-index proof.
"""

from fractions import Fraction
from math import isqrt

LIMIT = 2048
PRODUCT_LIMIT = 96
TAYLOR_LIMIT = 32

# Published data, independent of the algorithms below. A384819 starts at n=1.
SOURCE_A = (0, 1, 2, 1, 4, 3, 6, 1, 2, 7, 10, 3, 12, 11, 3, 1, 16, 3, 18, 15, 14, 19, 22, 3, 4, 23, 2, 27, 28, 12, 30, 1, 15, 31, 7, 3, 36, 35, 32, 7, 40, 37, 42, 7, 21, 43, 46, 3, 6, 7, 27, 11, 52, 3, 34, 35, 50, 55, 58, 44, 60, 59, 5, 1, 18, 0, 66, 19, 39, 40, 70, 3, 72, 71, 3, 23, 1, 19, 78, 55, 2, 79, 82, 41, 47, 83, 51, 47, 88, 84, 74, 31, 86, 91, 17, 3, 96, 11, 42, 15, 100, 18, 102, 67, 93, 103, 106, 3, 108, 105, 104, 67, 112, 43, 48, 43, 77, 115, 42, 12, 10, 119, 75, 47, 4, 64, 126, 1, 122, 37, 130, 24, 14, 131, 48, 107, 136, 36, 138, 28, 87, 139, 19, 3, 37, 143, 14, 59, 148, 12, 150, 127, 108, 46, 134, 67, 156, 155, 99, 87, 95, 3, 162, 67, 141, 163, 166, 1, 12, 60, 149, 71, 172, 54, 107, 31, 111, 175, 178, 152, 180, 75, 176, 167, 127, 91, 178, 79, 113, 116, 190, 3, 192, 191, 116, 27, 196, 9, 198, 7, 194, 199, 160, 114, 184, 203, 174, 83, 174, 78, 210, 91, 135, 211, 108, 3, 112, 215, 212, 121, 13, 115, 222, 195, 21, 223, 226, 163, 228, 137, 123, 227, 232, 82, 167, 103, 230, 121, 238, 204, 240, 19, 2, 107, 105, 90, 242, 247, 159, 7, 250, 212, 142, 251, 16, 1, 256, 139, 179, 57, 240, 259, 262, 208, 138, 5, 171, 119, 268, 57, 270, 187, 184, 271, 259, 204, 276, 275, 14, 164, 280, 108, 282, 127, 196, 172, 39, 3, 16, 226, 284, 131, 292, 37, 97, 11, 285, 295, 234, 44, 258, 299, 195, 239, 284, 189, 306, 178, 302, 155, 310, 51, 312, 311, 201, 143, 316, 126, 34, 215, 207, 160, 31, 3, 68, 323, 320, 19, 98, 14, 330, 151, 32, 331, 247, 177, 336, 23, 219, 40, 323, 268, 6, 23, 8, 343, 346, 294, 348, 90, 131, 319, 352, 144, 334, 163, 338, 355, 358, 264, 18, 359, 15, 23, 198, 211, 366, 343, 3, 240, 169, 355, 372, 283, 3, 31, 162, 118, 378, 52, 374, 379, 382, 3, 47, 383, 50, 179, 388, 165, 3, 35, 255, 391, 137, 357, 396, 395, 238, 55, 400, 235, 122, 187, 291, 163, 185, 226, 408, 235, 267, 191, 252, 369, 228, 147, 410, 219, 418, 210, 420, 419, 15, 43, 297, 180, 140, 199, 18, 337, 430, 3, 432, 373, 387, 203, 25, 259, 438, 233, 5, 214, 442, 7, 157, 443, 291, 259, 448, 84, 53, 211, 446, 451, 238, 339, 456, 455, 378, 397, 460, 264, 462, 35, 367, 463, 466, 202, 225, 330, 464, 55, 448, 283, 367, 365, 27, 475, 478, 332, 435, 479, 10, 7, 367, 3, 486, 59, 482, 432, 490, 474, 394, 411, 141, 55, 454, 216, 498, 15, 327, 499, 502, 244, 484, 501, 32, 239, 508, 149, 322, 1, 14, 511, 288, 31, 312, 370, 339, 201)
# Ordinary coefficients of exp(sum (n^2-a(n))*x^n/n), starting at n=0.
SOURCE_F = (1, 1, 2, 4, 8, 14, 25, 43, 74, 124, 205, 335, 543, 869, 1379, 2170, 3388, 5249, 8079, 12353, 18776, 28375, 42651, 63782, 94923, 140614, 207384, 304578, 445528, 649200, 942495, 1363447, 1965697, 2824676, 4046190, 5778273, 8227533, 11681632, 16540183, 23357053, 32898242)


def coefficient_construction(limit):
    """Use only F'=B'F and the allowed interval, without prime-power formulas."""
    a = [0] * (limit + 1)
    b = [0] * (limit + 1)
    f = [1] + [0] * limit
    for n in range(1, limit + 1):
        # n*f[n] = n^2 - a[n] + sum_{k<n} b[k]*f[n-k].
        numerator_before_correction = n*n + sum(b[k]*f[n-k] for k in range(1, n))
        a[n] = numerator_before_correction % n
        b[n] = n*n - a[n]
        numerator = numerator_before_correction - a[n]
        assert numerator % n == 0
        f[n] = numerator // n
        assert 0 <= a[n] < n
        if n <= 64:
            # Exhaust all allowed corrections to test finite uniqueness.
            valid = [c for c in range(n) if (numerator_before_correction-c) % n == 0]
            assert valid == [a[n]]
    return a, b, f


def ledger_construction(limit):
    """Independent triangular recursion using only proper-divisor ghost sums."""
    a = [0] * (limit + 1)
    exponents = [0] * (limit + 1)
    proper_sum = [0] * (limit + 1)
    for n in range(1, limit + 1):
        a[n] = (n*n - proper_sum[n]) % n
        numerator = n*n - a[n] - proper_sum[n]
        assert numerator % n == 0
        exponents[n] = numerator // n
        for multiple in range(2*n, limit + 1, n):
            proper_sum[multiple] += n * exponents[n]
    return a, exponents


def euler_product(exponents, degree):
    """Multiply (1-x^d)^(-e[d]) using exact generalized binomial coefficients."""
    product = [1] + [0] * degree
    for d in range(1, degree + 1):
        factor = [1]
        for j in range(degree // d):
            numerator = factor[-1] * (exponents[d] + j)
            assert numerator % (j + 1) == 0
            factor.append(numerator // (j + 1))
        old = product
        product = [sum(old[n-d*j] * factor[j] for j in range(n//d+1))
                   for n in range(degree + 1)]
    return product


def multiply(left, right, degree):
    return [sum(left[k]*right[n-k] for k in range(n+1))
            for n in range(degree+1)]


def taylor_exponential(b, degree):
    """Direct rational Taylor expansion, independent of the derivative recurrence."""
    exponent = [Fraction(0)] + [Fraction(b[n], n) for n in range(1, degree+1)]
    term = [Fraction(1)] + [Fraction(0)]*degree
    total = term.copy()
    for k in range(1, degree+1):
        term = [c/k for c in multiply(term, exponent, degree)]
        total = [x+y for x,y in zip(total, term)]
    return total


def prime(n):
    return n >= 2 and all(n%d for d in range(2, isqrt(n)+1))


def main():
    a, b, f = coefficient_construction(LIMIT)
    ledger_a, exponents = ledger_construction(LIMIT)
    assert a == ledger_a
    assert tuple(a[1:521]) == SOURCE_A
    assert tuple(f[:len(SOURCE_F)]) == SOURCE_F
    assert euler_product(exponents, PRODUCT_LIMIT) == f[:PRODUCT_LIMIT+1]
    assert taylor_exponential(b, TAYLOR_LIMIT) == f[:TAYLOR_LIMIT+1]

    # Necessary ghost identity, checked for every n in the finite range.
    for n in range(1, LIMIT+1):
        assert b[n] == sum(d*exponents[d] for d in range(1,n+1) if n%d == 0)

    # Negative control: replacing a(2)=1 by a(2)=0 gives [x^2]F=5/2.
    bad_b = [0, 1, 4]
    bad_f = taylor_exponential(bad_b, 2)
    assert bad_f[2] == Fraction(5,2) and bad_f[2].denominator != 1

    pairs = []
    for p in range(2, LIMIT+1):
        if not prime(p):
            continue
        power, previous, exponent = p, 1, 1
        while power <= LIMIT:
            assert a[power] == p-1, (p, exponent, a[power])
            assert b[power]-b[previous] == power*exponents[power]
            if exponent >= 2:
                assert (power*power-previous*previous) % power == 0
                assert a[power] == a[previous]
            pairs.append((p, exponent, power))
            previous, power, exponent = power, power*p, exponent+1

    print("FIDELITY A384819: all 520 published terms (n=1..520) match")
    print(f"FIDELITY A384820: all {len(SOURCE_F)} displayed coefficients (n=0..{len(SOURCE_F)-1}) match")
    print(f"INDEPENDENCE: coefficient recurrence and divisor-ledger recursion agree for n=1..{LIMIT}")
    print(f"EULER_PRODUCT: exact signed-binomial product agrees through degree {PRODUCT_LIMIT}")
    print(f"FORMAL_EXPONENTIAL: direct Fraction Taylor expansion agrees through degree {TAYLOR_LIMIT}")
    print(f"GHOST_IDENTITY: all n=1..{LIMIT} satisfy b[n]=sum_(d|n) d*e[d]")
    print("BOUNDED_UNIQUENESS: all allowed corrections checked for n=1..64")
    print("NEGATIVE_CONTROL: a(2)=0 gives nonintegral coefficient 5/2")
    print(f"PRIME_POWERS_SOURCE: {sum(n<=520 for p,r,n in pairs)} indices <=520 satisfy a(p^r)=p-1")
    print(f"PRIME_POWERS_EXTENDED: {len(pairs)} indices <={LIMIT} satisfy a(p^r)=p-1 and the ghost-difference identity")
    print("A384819_PREFIX:", ",".join(map(str,a[1:17])))
    print("A384820_PREFIX:", ",".join(map(str,f[:17])))
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
