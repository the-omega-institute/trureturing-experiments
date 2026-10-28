#!/usr/bin/env python3
"""Exact, offline checks for Peter Kagey's OEIS A329479 recurrence.

Source: https://oeis.org/A329479/internal, revision 25.
No recurrence is used to produce the polynomial counts.
"""
from itertools import product
from math import comb
import platform
import sympy as sp


SOURCE_TERMS = (
    0, 0, 0, 0, 1, 2, 6, 15, 30, 66, 121, 242, 462, 903, 1806,
    3570, 7225, 14450, 29070, 58311, 116622, 233586, 466489, 932978,
    1864590, 3727815, 7455630, 14908530, 29822521, 59645042,
    119301006, 238612935, 477225870, 954473586, 1908903481,
    3817806962, 7635526542,
)
SOURCE_SEVEN = {
    (7, 5, 3), (7, 6, 5, 4, 3, 2), (7, 5, 1), (7, 3, 1),
    (7, 6, 5, 4, 2, 1), (7, 6, 4, 3, 2, 1),
}


def eval_mod_three(coeffs, k):
    """Horner evaluation, with coefficients listed in ascending order."""
    value = 0
    for bit in reversed(coeffs):
        value = (k * value + bit) % 3
    return value


def enumerate_degree(n):
    """Enumerate every monic binary degree-n polynomial, including e_0."""
    return [bits + (1,) for bits in product((0, 1), repeat=n)
            if all(eval_mod_three(bits + (1,), k) == 0 for k in range(3))]


def evaluation_dp(limit):
    """Track the joint evaluations at 0, 1, 2, then force the leading bit."""
    states = {(0, 0, 0): 1}
    counts = []
    for n in range(limit + 1):
        increment = tuple(pow(k, n, 3) for k in range(3))
        counts.append(states.get(tuple((-v) % 3 for v in increment), 0))
        nxt = dict(states)
        for state, count in states.items():
            target = tuple((a + b) % 3 for a, b in zip(state, increment))
            nxt[target] = nxt.get(target, 0) + count
        states = nxt
    return counts


def main():
    print('Python', platform.python_version(), 'SymPy', sp.__version__)
    limit = 1000
    counts = evaluation_dp(limit)
    assert tuple(counts[1:38]) == SOURCE_TERMS
    seven = enumerate_degree(7)
    exponents = {tuple(i for i in range(7, -1, -1) if e[i]) for e in seven}
    assert exponents == SOURCE_SEVEN
    print('FIDELITY: all 37 printed source terms and all six degree-seven polynomials match')

    for n in range(1, 19):
        explicit = enumerate_degree(n)
        assert len(explicit) == counts[n], n
        for e in explicit:
            assert e[0] == 0
            assert sum(e[i] for i in range(1, n + 1) if i % 2 == 0) % 3 == 0
            assert sum(e[i] for i in range(1, n + 1) if i % 2 == 1) % 3 == 0
    print('ENUMERATION: independent Horner checks at degrees 1..18:', tuple(counts[1:19]))

    b = [sum(comb(m, j) for j in range(0, m + 1, 3))
         for m in range(limit // 2 + 2)]
    r = [b[n // 2] * b[(n + 1) // 2] for n in range(limit + 1)]
    for n in range(1, limit + 1):
        assert counts[n] == r[n] - r[n - 1], n
    correction = (2, 1, -1, -2, -1, 1)
    for m in range(len(b)):
        assert 3 * b[m] == 2 ** m + correction[m % 6], m
    print('BRIDGE: evaluation DP agrees with parity-binomial products at degrees 1..1000')
    print('BINOMIAL: 3*B_m = 2^m + c_m checked at m = 0..501')

    coefficients = (2, 2, -5, -2, 10, -4, -4, 8)
    for n in range(9, limit + 1):
        assert counts[n] == sum(c * counts[n - i]
                                for i, c in enumerate(coefficients, 1)), n
    print('RECURRENCE: exact polynomial counts satisfy the conjecture for 9 <= n <= 1000')

    t, z = sp.symbols('t z')
    cgf = sum(correction[m] * t ** m for m in range(6)) / (1 - t ** 6)
    assert sp.cancel(cgf - (2 - t) / (1 - t + t ** 2)) == 0
    sqgf = sum(correction[m] ** 2 * t ** m for m in range(6)) / (1 - t ** 6)
    prgf = sum(correction[m] * correction[(m + 1) % 6] * t ** m
               for m in range(6)) / (1 - t ** 6)
    assert sp.cancel(sqgf - 2 / (1 - t) - (2 + t) / (1 + t + t ** 2)) == 0
    assert sp.cancel(prgf - 1 / (1 - t) - (1 - t) / (1 + t + t ** 2)) == 0
    h9 = (1 / (1 - 4*t) + (4 - 4*t) / (1 - 2*t + 4*t**2)
          + 2 / (1 - t) + (2 + t) / (1 + t + t**2))
    j9 = (2 / (1 - 4*t) + (5 - 8*t) / (1 - 2*t + 4*t**2)
          + 1 / (1 - t) + (1 - t) / (1 + t + t**2))
    q = (1 + z) * (1 - 2*z) * (1 - z + z**2) * (1 - 2*z**2 + 4*z**4)
    q_expanded = 1 - 2*z - 2*z**2 + 5*z**3 + 2*z**4 - 10*z**5 + 4*z**6 + 4*z**7 - 8*z**8
    assert sp.expand(q - q_expanded) == 0
    generating = (1 - z) * (h9.subs(t, z**2) + z*j9.subs(t, z**2)) / 9 - 1
    assert sp.cancel(generating - (z**5 + 4*z**8) / q) == 0
    print('ALGEBRA: exact rational cancellation gives A(z) = (z^5 + 4*z^8)/Q(z)')
    print('Q(z) =', sp.expand(q))
    print('Finite tests and rational-function algebra do not constitute a Lean proof of the count bridge.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
