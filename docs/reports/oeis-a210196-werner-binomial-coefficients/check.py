#!/usr/bin/env python3
"""Exact finite checks of OEIS A210196 revision 21; this is not a proof."""

from math import comb


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def add(*polynomials):
    """Add integer polynomials stored as coefficient lists in ascending order."""
    result = [0] * max(map(len, polynomials))
    for polynomial in polynomials:
        for k, coefficient in enumerate(polynomial):
            result[k] += coefficient
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def twice_x(polynomial):
    return [0] + [2 * coefficient for coefficient in polynomial]


def coefficient(polynomial, k):
    return polynomial[k] if k < len(polynomial) else 0


def predicted(n, k):
    if k == 0:
        return 1
    return (2 ** (k + 1)) * (comb(n - 1, k) if k <= n - 1 else 0)


def main():
    # Literal source update, always using the two OLD polynomial values.
    u, v = [1], [1]
    pairs = []
    for n in range(1, 101):
        pairs.append((u, v))
        u, v = add(u, v, [1]), add(twice_x(u), twice_x(v), [1])

    # Fidelity: the five rows explicitly printed in OEIS, independently transcribed.
    printed = [[1], [1, 4], [1, 8, 8], [1, 12, 24, 16], [1, 16, 48, 64, 32]]
    require([v for _, v in pairs[:5]] == printed, 'source five-row fidelity failed')
    for n, row in enumerate(printed, 1):
        print(f'FIDELITY row n={n}: {row}')
    print('FIDELITY_OK: all five printed rows and the first three source polynomials')

    checked = 0
    for n, (_, v) in enumerate(pairs, 1):
        require(len(v) == n, f'unexpected support at n={n}')
        for k in range(n + 5):
            require(coefficient(v, k) == predicted(n, k), f'formula failed at n={n}, k={k}')
            checked += 1
        # An additional published result, checked directly on the recurrence output.
        alternating = sum(c if k % 2 == 0 else -c for k, c in enumerate(v))
        require(alternating == (1 if n % 2 == 1 else -3), f'alternating row sum at n={n}')
    print(f'IDENTITY_OK: n=1..100, k=0..n+4; {checked} exact integer comparisons')
    print('BOUNDARIES_OK: n=1, k=0, top coefficient, and five zero coefficients past each row')
    print('ALTERNATING_ROW_SUMS_OK: n=1..100 gives 1,-3,1,-3,...')
    print('SCOPE: finite integer calculation only; no inference to all n or all k')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
