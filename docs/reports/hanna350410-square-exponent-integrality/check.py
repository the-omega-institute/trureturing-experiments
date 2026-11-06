#!/usr/bin/env python3
"""Exact recurrence and divisibility checks for Hanna A350410--A350412."""

from fractions import Fraction
from math import factorial


OEIS = {
    "A350410": [-1, 1, 3, 12, 116, -534, -19083, 4653416,
                8687989872, -128647855586426, -17178731504646708274,
                22932060477876906089601744, 336748983792218322012229442319058,
                -59340140012235302782558826850918965528904,
                -135935979554634873389993183019614435570206186506348,
                4359616519227623971178976912185258573259921937219097427608884],
    "A350411": [-1, 1, 6, 76, 74, 26058, -7794632, -14388774920,
                213380711849598, 28487519401892735386,
                -38029127603565193657832846,
                -558442679082054075230555562505128,
                98405860538000983882824515617294391388622,
                225427456282611271896161167867656072505435638401710],
    "A350412": [1, 12, 830, 1867901, 251714851830,
                3696867014099083814, 8468768030682252554158546818,
                4074040763513889480730186336041282515610,
                533088609878310228401628064629048305794583364205468626,
                23856105278172150534783174192628043609641782797360374708050534816035626],
}


def mul(a, b, limit):
    out = [0] * (min(limit, len(a) + len(b) - 2) + 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            if i + j <= limit:
                out[i + j] += ai * bj
    return out


def power(a, exponent, limit):
    out = [1]
    base = a[: limit + 1]
    while exponent:
        if exponent & 1:
            out = mul(out, base, limit)
        exponent >>= 1
        if exponent:
            base = mul(base, base, limit)
    return out


def coeff_sum(series, cutoff, degree):
    total = Fraction(0)
    for m in range(cutoff + 1):
        p = power(series, m * m, degree)
        total += Fraction(p[degree] if degree < len(p) else 0, factorial(m))
    return total


def solve(cutoff):
    # b_0=-1, b_1=1.  The third cutoff is the same B=A-1 normalization.
    b = [-1, 1]
    divisibility = []
    for n in range(2, 14):
        m = cutoff(n)
        c = power(b + [0] * (n - len(b) + 1), m * m, n)[n]
        assert c % m == 0, (n, m, c)
        divisibility.append((n, m, c // m))
        raw = coeff_sum(b + [0] * (n - len(b) + 1), m, n)
        slope = Fraction((-1) ** (m - 1), factorial(m - 1))
        bn = -raw / slope
        assert bn.denominator == 1, (n, bn)
        b.append(bn.numerator)
    return b, divisibility


def main():
    cases = {
        "A350410": lambda n: n,
        "A350411": lambda n: n + 1,
        "A350412": lambda n: 2 * n,
    }
    for name, cutoff in cases.items():
        b, divisibility = solve(cutoff)
        # A350412 lists A's coefficients from n=1, while b=A-1 has b_0=-1.
        reported = b[1:] if name == "A350412" else b
        expected = OEIS[name]
        assert reported[: min(len(expected), 12)] == expected[: min(len(expected), 12)]
        assert len(reported) >= 12
        print(f"{name}: first 12 coefficients = {reported[:12]}")
        print(f"{name}: divisibility checks n=2..13 = {len(divisibility)}")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
