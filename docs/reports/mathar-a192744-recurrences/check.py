#!/usr/bin/env python3
"""Offline exact checks for the two Mathar conjectures in OEIS A192744.

Source: https://oeis.org/A192744, revision 16, 16 December 2015.
The source terms/examples are used under https://oeis.org/LICENSE.
Only Python 3 standard-library integer arithmetic is used.
"""

from math import factorial


SOURCE_TERMS = (
    1, 1, 3, 8, 29, 133, 762, 5215, 41257, 369032, 3676209,
    40333241, 483094250, 6271446691, 87705811341, 1314473334832,
    21017294666173, 357096406209005, 6424799978507178,
    122024623087820183, 2439706330834135361, 51219771117454755544,
)
LIMIT = 200


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def trim(coefficients):
    """Ascending integer coefficients; the empty list denotes zero."""
    result = list(coefficients)
    while result and result[-1] == 0:
        result.pop()
    return result


def remainder_monic(dividend, divisor):
    """Ordinary integer long division by an arbitrary monic polynomial."""
    result, divisor = trim(dividend), trim(divisor)
    require(bool(divisor) and divisor[-1] == 1, "monic divisor")
    while len(result) >= len(divisor):
        offset = len(result) - len(divisor)
        leading = result[-1]
        for j, coefficient in enumerate(divisor):
            result[offset + j] -= leading * coefficient
        result = trim(result)
    return result


def source_reduction(polynomial):
    """Literal A192232 iteration: p=x^2*t+r, then p=(x+1)*t+r."""
    result = trim(polynomial)
    while len(result) > 2:
        quotient, remainder = result[2:], result[:2]
        next_result = remainder + [0] * (len(result) - 2)
        for j, coefficient in enumerate(quotient):
            next_result[j] += coefficient
            next_result[j + 1] += coefficient
        result = trim(next_result)
    return result


def coordinates(remainder):
    return tuple((remainder + [0, 0])[:2])


def main():
    # Generate p_n by the source definition, without either conjectured recurrence.
    polynomial = [1]
    polynomials = []
    remainders = []
    for n in range(LIMIT + 1):
        if n:
            polynomial = [factorial(n)] + polynomial
        polynomials.append(polynomial[:])
        monic_remainder = remainder_monic(polynomial, [-1, -1, 1])
        literal_remainder = source_reduction(polynomial)
        require(monic_remainder == literal_remainder, f"reduction bridge n={n}")
        remainders.append(coordinates(monic_remainder))

    expected_polynomials = (
        [1], [1, 1], [2, 1, 1], [6, 2, 1, 1], [24, 6, 2, 1, 1]
    )
    expected_remainders = ((1, 0), (1, 1), (3, 2), (8, 5), (29, 13))
    require(tuple(polynomials[:5]) == expected_polynomials, "source polynomials")
    require(tuple(remainders[:5]) == expected_remainders, "source reductions")
    a = [constant for constant, _ in remainders]
    require(tuple(a[:len(SOURCE_TERMS)]) == SOURCE_TERMS, "22 source terms")
    print("FIDELITY: all 5 displayed source polynomials and reductions reproduced")
    print("FIDELITY: all 22 published A192744 terms reproduced (n=0..21)")
    print("a(0)..a(10):", ",".join(map(str, a[:11])))
    print(f"REDUCTION: literal substitution equals monic remainder at n=0..{LIMIT}")

    # A third construction uses the two integral coordinates of 1 and phi.
    pair = (1, 0)
    for n in range(LIMIT + 1):
        require(pair == remainders[n], f"independent coordinates n={n}")
        pair = (pair[1] + factorial(n + 1), pair[0] + pair[1])
    print(f"COORDINATES: independent pair construction agrees at n=0..{LIMIT}")

    g = {n: a[n] - a[n - 1] - a[n - 2] for n in range(2, LIMIT + 1)}
    for n, value in g.items():
        require(value == (n - 1) * factorial(n - 1), f"forcing n={n}")
    print(f"FORCING: g(n)=(n-1)*(n-1)! at n=2..{LIMIT} ({LIMIT - 1} cases)")

    for n in range(4, LIMIT + 1):
        residual = (a[n] - (n + 2) * a[n - 1] + 2 * (n - 1) * a[n - 2]
                    + 3 * a[n - 3] - (n - 2) * a[n - 4])
        require(residual == 0, f"2014 conjecture n={n}")
        require(g[n] - (n + 1) * g[n - 1] + (n - 2) * g[n - 2] == 0,
                f"2014 forcing elimination n={n}")
    print(f"2014 CONJECTURE: zero residual at every n=4..{LIMIT} ({LIMIT - 3} cases)")

    for n in range(3, LIMIT + 1):
        residual = ((2 - n) * a[n] + (n * n - n - 1) * a[n - 1]
                    + (-n * n + 3 * n - 3) * a[n - 2]
                    - (n - 1) ** 2 * a[n - 3])
        require(residual == 0, f"2015 conjecture n={n}")
        require((n - 2) * g[n] - (n - 1) ** 2 * g[n - 1] == 0,
                f"2015 forcing elimination n={n}")
    print(f"2015 CONJECTURE: zero residual at every n=3..{LIMIT} ({LIMIT - 2} cases)")
    print("ARITHMETIC: integers only; no floating point and no fitted recurrence")
    print("SCOPE: finite exact checks; no claim of a Lean proof or an unbounded proof")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
