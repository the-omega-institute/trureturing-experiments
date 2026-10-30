#!/usr/bin/env python3
"""Exact checks for Arizmendi--Arizmendi Conjecture 3.10 (arXiv:2012.14114v2).

No floating-point evaluations are used. Spectra come from actual adjacency
matrices; all Shapley coefficients are rational. Radical signs are certified
by integer square-root enclosures. The infinite-family certificate uses the
monotone integral inequalities explained in README.md, not finite sampling.
"""
from fractions import Fraction
from functools import lru_cache
from math import comb, factorial, isqrt
import sympy as s


def eq(a, b):
    assert s.simplify(a - b) == 0, (a, b)


def bipartite(r, k):
    return s.Matrix(r + k, r + k,
                    lambda i, j: int((i < r) != (j < r)))


@lru_cache(None)
def induced_energy(rows):
    A = s.Matrix(rows)
    if not rows:
        return s.S.Zero
    eigen = A.eigenvals()
    assert sum(eigen.values()) == A.rows
    assert all(x.is_real is True for x in eigen)
    return s.simplify(sum(s.Abs(x) * mult for x, mult in eigen.items()))


def all_worths(A):
    n = A.rows
    return {mask: induced_energy(tuple(tuple(A[i, j] for j in range(n) if mask >> j & 1)
                                     for i in range(n) if mask >> i & 1))
            for mask in range(1 << n)}


def direct_shapley(A):
    n = A.rows
    w = all_worths(A)
    phi = [s.simplify(sum(s.Rational(factorial(mask.bit_count()) *
                                   factorial(n - mask.bit_count() - 1), factorial(n)) *
                          (w[mask | (1 << i)] - w[mask])
                          for mask in range(1 << n) if not mask >> i & 1))
           for i in range(n)]
    eq(sum(phi), w[(1 << n) - 1])
    return w, phi


def interval(expr, digits=16):
    """Closed rational enclosure of a rational linear combination of sqrt(N)."""
    scale = 10 ** digits
    lo = hi = Fraction(0)
    for term in s.Add.make_args(s.expand(expr)):
        coeff, root = term.as_coeff_Mul()
        assert coeff.is_Rational
        c = Fraction(int(s.numer(coeff)), int(s.denom(coeff)))
        if root == 1:
            a = b = Fraction(1)
        else:
            assert root.is_Pow and root.exp == s.Rational(1, 2)
            rad = int(root.base)
            assert rad >= 0 and root.base == rad
            z = isqrt(rad * scale * scale)
            assert z * z <= rad * scale * scale < (z + 1) * (z + 1)
            a = Fraction(z, scale)
            b = a if z * z == rad * scale * scale else Fraction(z + 1, scale)
        lo += c * (a if c >= 0 else b)
        hi += c * (b if c >= 0 else a)
    assert lo <= hi
    return lo, hi


def phi_a(m):
    H = (m + 1) * (m + 2)
    q = s.sqrt(2)
    return s.expand(2 * sum((m + q + (q - 2) * k) * s.sqrt(k)
                             for k in range(m + 1)) / H)


def phi_large_part(m):
    """Independent coalition-count formula for one B vertex, not efficiency."""
    n = m + 2
    total = s.S.Zero
    for r in range(3):
        for k in range(m):
            coefficient = s.Rational(comb(2, r) * comb(m - 1, k) *
                                     factorial(r + k) * factorial(n - r - k - 1),
                                     factorial(n))
            total += coefficient * 2 * (s.sqrt(r * (k + 1)) - s.sqrt(r * k))
    return s.expand(total)


# Fidelity: Example 3.9, with the center indexed 0 instead of 2.
w, phi = direct_shapley(bipartite(1, 2))
expected = [(2 + 2 * s.sqrt(2)) / 3] + [(2 * s.sqrt(2) - 1) / 3] * 2
for x, y in zip(phi, expected):
    eq(x, y)
for mask, value in w.items():
    slack = sum(phi[i] for i in range(3) if mask >> i & 1) - value
    assert interval(slack)[0] >= 0
print('FIDELITY Example 3.9: center=(2+2*sqrt(2))/3; leaves=(2*sqrt(2)-1)/3; all 8 core constraints pass')

# Actual adjacency characteristic polynomials, including all empty-part cases.
z = s.Symbol('z')
for r in range(5):
    for k in range(6):
        A = bipartite(r, k)
        expected_poly = z ** (r + k) if r * k == 0 else z ** (r + k - 2) * (z*z - r*k)
        assert s.expand(A.charpoly(z).as_expr() - expected_poly) == 0
print('SPECTRAL: 30 exact characteristic polynomials of K_(r,k), 0<=r<=4, 0<=k<=5')

# Independent literal subset summation versus the reduced family expression.
for m in range(1, 5):
    _, ph = direct_shapley(bipartite(2, m))
    eq(ph[0], phi_a(m))
    eq(ph[1], phi_a(m))
    for p in ph[2:]:
        eq(p, phi_large_part(m))
print('COUNTING: full subset Shapley sums agree with both orbit formulas for K_(2,m), 1<=m<=4')

# Exact factorial-to-probability bridge for both predecessor classes.
for m in range(71):
    H = (m + 1) * (m + 2)
    mass = Fraction(0)
    for k in range(m + 1):
        p0 = Fraction(comb(m, k) * factorial(k) * factorial(m + 1 - k), factorial(m + 2))
        p1 = Fraction(comb(m, k) * factorial(k + 1) * factorial(m - k), factorial(m + 2))
        assert p0 == Fraction(m + 1 - k, H)
        assert p1 == Fraction(k + 1, H)
        mass += p0 + p1
    assert mass == 1
print('PROBABILITIES: both factorial identities and total mass 1 pass for 0<=m<=70')

# Infinite-family certificate: exact algebra after monotone integral comparison.
m, q, t, x = s.symbols('m q t x', positive=True)
H = (m + 1) * (m + 2)
# Integrals are exact, with their monotone comparison directions kept explicit.
eq(s.integrate(s.sqrt(x), (x, 0, m)), s.Rational(2, 3) * m ** s.Rational(3, 2))
eq(s.integrate(x ** s.Rational(3, 2), (x, 0, m)), s.Rational(2, 5) * m ** s.Rational(5, 2))
P = (13 - 9*q)*m*m/15 + (1 - 4*q/3)*m + 2 - 2*q
bound = ((m + q)*s.Rational(2, 3)*m + (q - 2)*(s.Rational(2, 5)*m*m + m)) / H - (q - 1)
eq(bound, P/H)
assert Fraction(99, 70) ** 2 - 2 == Fraction(1, 4900)
assert Fraction(99, 70) < 2
R = (19*m*m - 930*m - 870)/1050
eq(P.subs(q, s.Rational(99, 70)), R)
eq(P - R, (s.Rational(99, 70) - q) * (s.Rational(3, 5)*m*m + s.Rational(4, 3)*m + 2))
shifted = s.Poly(s.expand(1050 * R.subs(m, t + 50)), t)
assert shifted.all_coeffs() == [19, 970, 130]
assert all(c > 0 for c in shifted.all_coeffs())
print('FAMILY CERTIFICATE: q=sqrt(2)<99/70, since (99/70)^2-2=1/4900')
print('FAMILY CERTIFICATE: Delta_m/(2*sqrt(m)) >= P(m)/((m+1)*(m+2)) > (19*m^2-930*m-870)/(1050*(m+1)*(m+2))')
print('FAMILY CERTIFICATE: m=t+50 gives 19*t^2+970*t+130>0 for every t>=0')

# Direct, rigorously enclosed deficits; independently calculated B payoff.
for mm in [50, 64, 100]:
    pa, pb = phi_a(mm), phi_large_part(mm)
    eq(2*pa + mm*pb, 2*s.sqrt(2*mm))
    delta = pa - 2*(s.sqrt(2)-1)*s.sqrt(mm)
    coalition_deficit = 2*s.sqrt(mm) - (pa + mm*pb)
    eq(delta, coalition_deficit)
    low, high = interval(delta)
    assert low > 0
    print(f'COUNTEREXAMPLE K_(2,{mm}): exact deficit enclosure [{low}, {high}], lower endpoint > 0')

# Explicit 52-vertex graph certificate independent of the energy shortcut.
A = bipartite(2, 50)
assert A == A.T and all(A[i, i] == 0 for i in range(52))
assert set(A) == {s.S.Zero, s.S.One}
assert min(sum(A[i, j] for j in range(52)) for i in range(52)) == 2
seen, frontier = {0}, [0]
while frontier:
    i = frontier.pop()
    for j in range(52):
        if A[i, j] and j not in seen:
            seen.add(j)
            frontier.append(j)
assert len(seen) == 52
assert A.rank() == 2
assert A*A*A == 100*A and s.trace(A*A) == 200
print('ACTUAL GRAPH K_(2,50): 52 vertices, simple, connected (complete bipartite with nonempty parts), minimum degree 2, rank 2, A^3=100*A, tr(A^2)=200')

# Asymptotic coefficient from the two monotone Riemann-sum limits.
limit = s.simplify(2*(s.Rational(2, 3) + (s.sqrt(2)-2)*s.Rational(2, 5)) - 2*(s.sqrt(2)-1))
eq(limit, (26 - 18*s.sqrt(2))/15)
assert interval(limit)[0] > 0
assert 26*26 - 18*18*2 == 28
print('ASYMPTOTIC: Delta_m/sqrt(m) -> (26-18*sqrt(2))/15 > 0; exact positivity: 26^2-18^2*2=28')
print('ALL_CHECKS_OK')
