#!/usr/bin/env python3
"""Exact, finite fidelity and arithmetic checks for OEIS A227051.

Source: https://oeis.org/A227051, revision 10.
This finite check is not a proof of the conjecture for unbounded n.
Python 3 standard library only; no network access or floating point.
"""
from fractions import Fraction as Q
from math import comb, factorial

LIMIT = 120
SOURCE_DATA = (
    1, 3, 53, 2333, 191673, 25307913, 4900979093, 1308599657693,
    460755584003313, 206844794964734673, 115313659955341400333,
    78158334287649490486853, 63294640267864707577746153,
    60357724113527363258814802233, 66943314938342593826952764256773,
    85443499990582824984241143043808813,
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    degree = 2 * LIMIT
    facts = [factorial(k) for k in range(degree + 1)]
    # Literal ordinary rational Taylor series, independently of the proposed
    # closed denominator formula and of the integer inverse recurrence.
    cosine = [Q((-1) ** (k // 2), facts[k]) if k % 2 == 0 else Q(0)
              for k in range(degree + 1)]
    sine = [Q((-1) ** (k // 2), facts[k]) if k % 2 else Q(0)
            for k in range(degree + 1)]
    sinh = [Q(1, facts[k]) if k % 2 else Q(0)
            for k in range(degree + 1)]
    denominator = [cosine[k] - sum((sine[j] * sinh[k - j]
                    for j in range(k + 1)), Q(0))
                   for k in range(degree + 1)]
    require(denominator[0] == 1, "denominator constant term")
    inverse = [Q(1)]
    for k in range(1, degree + 1):
        inverse.append(-sum((denominator[j] * inverse[k - j]
                              for j in range(1, k + 1)), Q(0)))
    for k in range(1, degree + 1):
        require(sum((denominator[j] * inverse[k - j]
                     for j in range(k + 1)), Q(0)) == 0,
                f"ordinary inverse product at degree {k}")
    require(all(inverse[k] == denominator[k] == 0
                for k in range(1, degree + 1, 2)), "odd coefficients")
    normalized = [facts[2*n] * inverse[2*n] for n in range(LIMIT + 1)]
    require(tuple(normalized[:len(SOURCE_DATA)]) == SOURCE_DATA,
            "all 16 published OEIS DATA values")
    print(f"FIDELITY: reproduced all {len(SOURCE_DATA)} OEIS DATA values exactly")
    print("INITIAL: " + ",".join(str(x) for x in SOURCE_DATA[:7]))
    print(f"FORMAL_SERIES: product and odd-zero checks through degree {degree}")

    # Separately derive factorial-normalized denominator coefficients by
    # finite binomial sums and invert by an integer triangular recurrence.
    d = [(-1)**n - sum(comb(2*n, 2*k+1) * (-1)**k for k in range(n))
         for n in range(LIMIT + 1)]
    a = [1]
    for n in range(1, LIMIT + 1):
        a.append(-sum(comb(2*n, 2*j) * d[j] * a[n-j]
                      for j in range(1, n + 1)))
    for n in range(LIMIT + 1):
        require(d[n] == facts[2*n] * denominator[2*n], f"denominator bridge {n}")
        expected_d = 1 if n % 2 == 0 else -1 - 2 * (-4)**((n-1)//2)
        require(d[n] == expected_d, f"closed denominator formula {n}")
        require(normalized[n].denominator == 1, f"integrality {n}")
        require(a[n] == normalized[n], f"independent inverse bridge {n}")
        require(a[n] % 2 == 1, f"mod 2 {n}")
        if n == 0:
            require(a[n] == 1 and a[n] % 10 != 3, "excluded index zero")
            continue
        require(a[n] % 5 == 3 and a[n] % 10 == 3, f"last digit {n}")
        t = sum(comb(2*n, 2*j) for j in range(n+1))
        e = sum((-1)**j * comb(2*n, 2*j) for j in range(n+1))
        s = sum(comb(2*n, 2*j) * d[j] for j in range(n+1))
        require(t == 2**(2*n-1), f"even binomial sum {n}")
        require(e == (0 if n % 2 else (-4)**(n//2)), f"alternating sum {n}")
        require((t % 5, e % 5, s % 5) ==
                ((2, 0, 3) if n % 2 else (3, 1, 4)), f"mod 5 sums {n}")
        require((3*s - 2*d[n]) % 5 == 0, f"inverse mod 5 cancellation {n}")
        require(sum(comb(2*n, 2*j) * (1 if n-j == 0 else 3) * d[j]
                    for j in range(n+1)) % 5 == 0, f"candidate inverse mod 5 {n}")
    print(f"BRIDGES: rational Taylor inversion = integer binomial recurrence for 0 <= n <= {LIMIT}")
    print(f"DENOMINATOR: d(2r)=1; d(2r+1)=-1-2*(-4)^r checked for 0 <= n <= {LIMIT}")
    print(f"RESIDUES: a(n) odd for 0 <= n <= {LIMIT}; a(n) == 3 mod 5 and mod 10 for 1 <= n <= {LIMIT}")
    print(f"IDENTITIES: even/alternating binomial sums and mod-5 inverse cancellation for 1 <= n <= {LIMIT}")
    # w_6 = 6 [x^6] log(A) = [x^5] A'/A = -[x^5] D'/D.
    w6 = -sum((k * denominator[k] * inverse[6-k] for k in range(1, 7)), Q(0))
    require(w6 == Q(379, 60), "logarithmic weight screen")
    print(f"LITERATURE_HYPOTHESIS: w_6 = {w6}; denominator divisible by both 2 and 5")
    print("BOUNDARY: a(0)=1, deliberately excluded from the conjecture")
    print("ALL_CHECKS_OK")


if __name__ == '__main__':
    main()
