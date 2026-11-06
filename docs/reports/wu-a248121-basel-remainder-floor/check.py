#!/usr/bin/env python3
"""Exact checks for OEIS A248121; finite evidence, not a Lean proof.

Source: https://oeis.org/A248121, revision 18, NAME, DATA and FORMULA.
The 48 DATA terms are reproduced below with OEIS attribution.
Arithmetic uses Python integers, Fraction and exact SymPy expressions only.
"""
from fractions import Fraction as F
import sympy as S

SOURCE_TERMS = (
    2, 9, 20, 34, 53, 76, 102, 133, 168, 206, 249, 296, 346, 401,
    460, 522, 589, 660, 734, 813, 896, 982, 1073, 1168, 1266, 1369,
    1476, 1586, 1701, 1820, 1942, 2069, 2200, 2334, 2473, 2616,
    2762, 2913, 3068, 3226, 3389, 3556, 3726, 3901, 4080, 4262, 4449, 4640,
)
LIMIT = 1000


def atan_reciprocal_bounds(q, terms=64):
    """Alternating-series rational enclosure for atan(1/q), q > 1."""
    assert q > 1 and terms > 0
    partial = sum((F((-1) ** k, (2 * k + 1) * q ** (2 * k + 1))
                   for k in range(terms)), F(0))
    following = partial + F((-1) ** terms, (2 * terms + 1) * q ** (2 * terms + 1))
    return min(partial, following), max(partial, following)


def predicted(n):
    return 2 * n * n + (2 * n) // 3


def main():
    # Machin's classical identity: pi = 16 atan(1/5) - 4 atan(1/239).
    a5lo, a5hi = atan_reciprocal_bounds(5)
    a239lo, a239hi = atan_reciprocal_bounds(239)
    pi_lo = 16 * a5lo - 4 * a239hi
    pi_hi = 16 * a5hi - 4 * a239lo
    assert 3 < pi_lo < pi_hi < F(22, 7)
    zeta_lo, zeta_hi = pi_lo * pi_lo / 6, pi_hi * pi_hi / 6
    harmonic2 = F(0)
    values = [0]
    for n in range(1, LIMIT + 1):
        harmonic2 += F(1, n * n)
        d_lo = F(1, n) + harmonic2 - zeta_hi
        d_hi = F(1, n) + harmonic2 - zeta_lo
        assert 0 < d_lo < d_hi
        lo, hi = 1 / d_hi, 1 / d_lo
        expected = predicted(n)
        # The entire certified reciprocal interval lies in one floor cell.
        assert expected <= lo < hi < expected + 1, n
        values.append(expected)
        if n <= len(SOURCE_TERMS):
            assert expected == SOURCE_TERMS[n - 1], n
    print('FIDELITY: all 48 published OEIS DATA terms reproduced from certified rational pi bounds')
    print('SETTLEMENT: original analytic definition has floor 2*n^2+floor(2*n/3) for n=1..1000')
    print('FIRST_TERMS:', ','.join(map(str, values[1:13])))

    n = S.symbols('n', positive=True, integer=True)
    h = 1 / (n * (n + 1) ** 2)
    upper = 3 / (2 * n * (3 * n + 1))
    lower = 9 / (2 * (9 * n ** 2 + 3 * n + 1))
    upper_gap = 2 / (n * (n + 1) ** 2 * (3 * n + 1) * (3 * n + 4))
    lower_gap = -(6 * n + 13) / (
        n * (n + 1) ** 2 * (9 * n ** 2 + 3 * n + 1) * (9 * n ** 2 + 21 * n + 13))
    assert S.cancel(h - (1/n - 1/(n+1) - 1/(n+1)**2)) == 0
    assert S.cancel(upper - upper.subs(n, n+1) - h - upper_gap) == 0
    assert S.cancel(lower - lower.subs(n, n+1) - h - lower_gap) == 0
    base = 2*n**2 + S.Rational(2, 3)*n
    assert S.cancel(1/upper - base) == 0
    assert S.cancel(1/lower - base - S.Rational(2, 9)) == 0
    for residue in range(3):
        fractional = F(2 * residue, 3) % 1
        assert fractional + F(2, 9) < 1
    print('SYMBOLIC: tail summand identity, both signed barrier identities and inverse width 2/9 verified')
    print('RESIDUES: fractional parts 0,2/3,1/3 each remain below 1 after adding 2/9')

    for n in range(6, LIMIT + 1):
        assert values[n] == (2*values[n-1] - values[n-2] + values[n-3]
                             - 2*values[n-4] + values[n-5])
    x = S.symbols('x')
    conjectured = x*(1+x)**2*(x+2)/((1-x)**3*(1+x+x*x))
    derived = 2*x*(1+x)/(1-x)**3 + (x*x+x**3)/((1-x)*(1-x**3))
    assert S.cancel(conjectured - derived) == 0
    denominator = [1, -2, 1, -1, 2, -1]
    numerator = [0, 2, 5, 4, 1]
    assert S.expand((1-x)**3*(1+x+x*x)) == sum(c*x**k for k,c in enumerate(denominator))
    for n in range(LIMIT + 1):
        convolution = sum(denominator[j]*values[n-j] for j in range(min(n, 5)+1))
        assert convolution == (numerator[n] if n < len(numerator) else 0), n
    print('RECURRENCE: all n=6..1000 verified; generating-function coefficient identity through degree 1000')
    print('GENERATING_FUNCTION: exact rational-function simplification verified')
    print('SCOPE: checks certify finite cases and algebraic identities; infinite-tail bounds still require proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
