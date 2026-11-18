#!/usr/bin/env python3
"""Exact copula-CDF integration and a separating certificate; no numerical evidence."""
import sympy as s

R = s.Rational
y, a, q, t, n = s.symbols('y a q t n', real=True)

def equal(actual, expected):
    assert s.factor(actual - expected) == 0, (actual, expected)

def integrate_law(xwidths, ywidths, density):
    """Piecewise-constant joint density, ordered positive widths summing to one."""
    equal(sum(xwidths), 1)
    equal(sum(ywidths), 1)
    for row in density:
        equal(sum(c * w for c, w in zip(row, ywidths)), 1)
    for j in range(len(ywidths)):
        equal(sum(xwidths[i] * density[i][j] for i in range(len(xwidths))), 1)
    cdf_square = 0
    mean_square = 0
    copy_distance = 0
    means = []
    for wx, row in zip(xwidths, density):
        left = prefix = mean = row_cdf_square = row_distance = 0
        for wy, c in zip(ywidths, row):
            right = left + wy
            F = prefix + c * (y - left)
            row_cdf_square += s.integrate(F**2, (y, left, right))
            row_distance += 2 * s.integrate(F * (1 - F), (y, left, right))
            mean += c * (right**2 - left**2) / 2
            prefix += c * wy
            left = right
        equal(prefix, 1)
        means.append(s.factor(mean))
        cdf_square += wx * row_cdf_square
        mean_square += wx * mean**2
        copy_distance += wx * row_distance
    xi = s.factor(6 * cdf_square - 2)
    eta = s.factor(12 * (mean_square - R(1, 4)))
    equal(xi, 1 - 3 * copy_distance)
    return xi, eta, means

# Source fidelity: the actual binary conditional CDFs G_sigma(v)=v+sigma min(v,a,1-v).
# The three response-cell widths are nonnegative for 0 <= a <= 1/2;
# zero-width boundary cells contribute zero integrals.
xi0, eta0, _ = integrate_law([R(1, 2)] * 2, [a, 1 - 2*a, a],
                            [[0, 1, 2], [2, 1, 0]])
b = 2*a**2*(3-4*a)
d = 12*a**2*(1-a)**2
equal(xi0, b)
equal(eta0, d)
print('SOURCE_FIDELITY binary model: xi =', xi0, '; eta =', eta0)

# Independent computation from the full joint density, with symbolic cell widths.
xw = [1-q, q/2, q/2]
yw = [1-q, q/3, q/3, q/3]
mat = [[1/(1-q), 0, 0, 0], [0, 2/q, 1/q, 0], [0, 0, 1/q, 2/q]]
xi, eta, means = integrate_law(xw, yw, mat)
equal(xi, 2*q-R(44, 27)*q**2)
equal(eta, 3*q-3*q**2+R(16, 27)*q**3)
print('FAMILY uniform marginals verified algebraically for 0 < q < 1')
print('FAMILY xi =', s.expand(xi), '; eta =', s.expand(eta))
q0 = R(3, 5)
equal(xi.subs(q, q0), R(46, 75))
equal(eta.subs(q, q0), R(106, 125))
assert [v.subs(q, q0) for v in means] == [R(1, 5), R(17, 30), R(5, 6)]
print('WITNESS q = 3/5: xi = 46/75; eta = 106/125; means = 1/5, 17/30, 5/6')

# Uniform separation of every seed branch, certified by polynomial identities
# and the displayed nonnegative factors (no finite-n extrapolation).
slope, M = R(4, 5), R(891, 2500)
X = 1-(1-b)/n
Y = 1-(1-d)/n**2
L = Y-slope*X
bracket1 = R(7, 20)+R(7, 5)*a-3*a**2
bracket2 = -3*a**2+R(8, 5)*a+R(4, 25)
equal(R(7, 20)-L.subs(n, 1), (1-2*a)**2 * bracket1)
equal(M-L.subs(n, 2), (a-R(1, 5))**2 * bracket2)
equal(bracket1, R(7, 20)-a/10+3*a*(R(1, 2)-a))
equal(bracket2, R(4, 25)+a/10+3*a*(R(1, 2)-a))
assert R(7, 20)-R(1, 2)/10 == R(3, 10) > 0
assert R(4, 25) > 0
equal(2*b-d, 4*a**3*(2-3*a))
base = R(1, 5)+R(4, 5)/n-1/n**2
equal(base-L, b*(4*n-10)/(5*n**2)+(2*b-d)/n**2)
equal(R(16, 45)-base, (R(1, 3)-1/n)*(R(7, 15)-1/n))
# On 0<=a<=1/2: b>=0 and 2b-d>=0. For n>=3 each preceding factor is >=0.
assert 3-4*R(1, 2) > 0 and 2-3*R(1, 2) > 0
assert 4*3-10 > 0 and R(7, 15)-R(1, 3) > 0
assert M-R(7, 20) == R(4, 625) > 0
assert M-R(16, 45) == R(19, 22500) > 0
assert M-(1-slope) > 0
print('ALL_BRANCHES separator y <= (4/5)x + 891/2500: n=1, n=2, every n>=3, and (1,1)')
# The halfspace is convex, so the same inequality holds on conv(S_u).
gap = s.factor(eta-slope*xi-M)
equal(gap.subs(q, q0), R(7, 7500))
equal(gap.subs(q, t+q0), R(7, 7500)+t/225-R(17, 27)*t**2+R(16, 27)*t**3)
lower = R(7, 7500)-R(1, 22500)-R(17, 270000)-R(16, 27000000)
assert lower == R(619, 750000) == R(1857, 2250000) > R(1, 1250) > 0
assert R(59, 100) > 0 and R(61, 100) < 1
print('INTERVAL 59/100 <= q <= 61/100: separation >=', lower, '> 1/1250')

# Everywhere-positive density: recompute the conditional CDFs of the mixture.
r = R(999, 1000)
positive_mat = [[r*c+1-r for c in row] for row in mat]
xip, etap, _ = integrate_law(xw, yw, positive_mat)
equal(xip, r**2*xi)
equal(etap, r**2*eta)
positive_gap = (etap-slope*xip-M).subs(q, q0)
assert positive_gap == R(13689, 62500000) > 0
print('POSITIVE_DENSITY c >= 1/1000: separation =', positive_gap)
print('ALL_CHECKS_OK')
