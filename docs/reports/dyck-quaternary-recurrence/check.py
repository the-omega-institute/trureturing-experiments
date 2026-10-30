#!/usr/bin/env python3
"""Exact checks for Mathar's OEIS A379464 recurrence; no Lean proof.

Source: https://oeis.org/A379464/internal, revision 20.
Python 3 and SymPy are required. All evidence uses integers, Fraction, or
symbolic polynomials over ZZ. Bounds are explicit and are not a proof that
all path counts satisfy the recurrence.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import comb, factorial
import sys

import sympy as sp

MAX_N = 512
SOURCE_TERMS = [
    1, 1, 1, 1, 4, 16, 46, 106, 226, 514, 1306, 3466, 9002, 22634,
    56330, 142026, 364743, 945303, 2448511, 6323695, 16336885,
    42363693, 110340297, 288229377, 753920796, 1973799396, 5174280216,
    13588243696, 35748326836, 94188788164, 248464963876, 656148369796,
]


def path_counts(bound, steps):
    """Count every admissible word by independent height-state DP."""
    states = {0: 1}
    result = [1]
    for _ in range(bound):
        nxt = defaultdict(int)
        for height, count in states.items():
            for step in steps:
                if height + step >= -2:
                    nxt[height + step] += count
        states = nxt
        result.append(states.get(0, 0))
    return result


def is_path(word):
    height = 0
    for step in word:
        height += step
        if height < -2:
            return False
    return height == 0


def b(k):
    value = Fraction(3 * comb(4 * k + 3, k), 4 * k + 3)
    assert value.denominator == 1
    assert value == Fraction(3 * factorial(4 * k + 2), factorial(k) * factorial(3 * k + 3))
    return value.numerator


def coefficients(n):
    """Transcribed, factored OEIS coefficients in their displayed order."""
    return [
        3 * (n + 4) * (3 * n + 4) * (3 * n + 8),
        3 * (-63 * n**3 - 297 * n**2 - 349 * n - 60),
        3 * (189 * n**3 + 270 * n**2 - 229 * n - 140),
        15 * (-63 * n**3 + 117 * n**2 + 44 * n - 64),
        689 * n**3 - 5372 * n**2 + 6946 * n - 1288,
        (n - 4) * (201 * n**2 + 2767 * n - 3011),
        -(n - 5) * (579 * n + 257) * (n - 4),
        229 * (n - 5) * (n - 6) * (n - 4),
    ]


def e_action(poly, m):
    """E(z^(m+j)) = (m+j)z^(m+j) - (m+j+1)z^(m+j+1)."""
    result = defaultdict(lambda: sp.Integer(0))
    for j, value in poly.items():
        result[j] += (m + j) * value
        result[j + 1] -= (m + j + 1) * value
    return dict(result)


def symbolic_certificate():
    m, n, x = sp.symbols('m n x', integer=True)
    powers = [{0: sp.Integer(1)}]
    for _ in range(3):
        powers.append(e_action(powers[-1], m))
    result = defaultdict(lambda: sp.Integer(0))
    # Left multiplication by (1-z)^4; no commuting with E.
    for power, scalar in enumerate([384, 528, 216, 27]):
        for degree, value in powers[power].items():
            for shift in range(5):
                result[degree + shift] += scalar * value * (-1)**shift * comb(4, shift)
    # Left multiplication by -256 z^4.
    for power, scalar in enumerate([90, 63, 14, 1]):
        for degree, value in powers[power].items():
            result[degree + 4] -= 256 * scalar * value
    assert set(result) == set(range(8))
    for j, target in enumerate(coefficients(n)):
        assert sp.Poly(sp.expand(result[j].subs(m, n-j) - target), n, domain=sp.ZZ).is_zero
    # The two cubics in the theta-to-E conjugation, multiplied by 64.
    assert sp.expand(64 * 3 * (x/4+1) * (3*x/4+1) * (3*x/4+2)
                     - (27*x**3+216*x**2+528*x+384)) == 0
    assert sp.expand(64 * 8 * (2*x/4+3) * (4*x/4+3) * (4*x/4+5)
                     - 256 * (x**3+14*x**2+63*x+90)) == 0
    return result


def check_forest_encoding():
    # Independent enumeration of all candidate preorder forest words at k <= 4.
    total_forests = 0
    for k in range(5):
        forests = set()
        for internal_positions in combinations(range(4*k+3), k):
            internal = set(internal_positions)
            word = tuple(3 if i in internal else -1 for i in range(4*k+3))
            slots = 3
            valid = True
            for i, step in enumerate(word):
                if slots <= 0:
                    valid = False
                    break
                slots += step
                if i < len(word)-1 and slots <= 0:
                    valid = False
            if valid and slots == 0:
                assert word[-3:] == (-1, -1, -1)
                path = tuple(-s for s in reversed(word[:-3]))
                assert is_path(path)
                assert path.count(-3) == k
                assert tuple(-s for s in reversed(path)) + (-1, -1, -1) == word
                forests.add(path)
        assert len(forests) == b(k)
        total_forests += len(forests)
    return total_forests


def main():
    counts = path_counts(MAX_N, (1, 0, -3))
    assert counts[:len(SOURCE_TERMS)] == SOURCE_TERMS
    paths4 = {tuple(w) for w in product((1, 0, -3), repeat=4) if is_path(w)}
    assert paths4 == {(0, 0, 0, 0), (1, -3, 1, 1), (1, 1, -3, 1), (1, 1, 1, -3)}
    print('fidelity: OEIS revision 20 terms n=0..31 and all four source n=4 paths agree')
    for n in range(11):
        assert sum(is_path(w) for w in product((1, 0, -3), repeat=n)) == counts[n]
    print('independent brute-force word enumeration: n=0..10 agrees with height DP')
    for n in range(MAX_N + 1):
        assert counts[n] == sum(comb(n, 4*k) * b(k) for k in range(n//4 + 1))
    nonzero = path_counts(MAX_N, (1, -3))
    for length, value in enumerate(nonzero):
        assert value == (b(length//4) if length % 4 == 0 else 0)
    print(f'exact path DP, zero insertion, and Scheuerle source sum agree: n=0..{MAX_N}')
    print(f'preorder forest encoding exhaustively checked: k=0..4; {check_forest_encoding()} forests')
    for k in range(1, MAX_N + 1):
        assert 3*(k+1)*(3*k+1)*(3*k+2)*b(k) == 8*(2*k+1)*(4*k-1)*(4*k+1)*b(k-1)
    assert 3 * 1 * 1 * 2 * b(0) == 6
    print(f'factorial-ratio recurrence: k=1..{MAX_N}; inhomogeneous theta constant = 6')
    symbolic_certificate()
    print('symbolic ZZ polynomial certificate: all eight coefficient polynomials agree identically')
    for n in range(MAX_N + 1):
        residual = sum(c*counts[n-j] for j, c in enumerate(coefficients(n)) if j <= n)
        forcing = 384*(-1)**n*comb(3,n) if n <= 3 else 0
        assert residual == forcing
    print(f'inhomogeneous L A = 384(1-z)^3: coefficients n=0..{MAX_N} agree')
    print(f'Mathar recurrence residual = 0: all {MAX_N-6} indices n=7..{MAX_N}')
    print('scope: finite path checks and exact operator algebra; the all-n counting bridge is not kernel-proved')
    print(f'Python {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}; SymPy {sp.__version__}')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
