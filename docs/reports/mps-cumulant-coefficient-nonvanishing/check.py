#!/usr/bin/env python3
"""Exact checks for MPS arXiv:2311.04831v2, Remark 3.31 (no Lean proof).

Requires Python 3 and SymPy. All certificates use integers and rational algebra.
The source's R recurrence independently generates the needed coefficient forms.
"""
import sympy as s

p, t, sigma, b = s.symbols('p t sigma b', real=True)
T = s.symbols('T0:15')
Q = s.Rational


def same(x, y=0):
    assert s.expand(x - y) == 0, (x, y)


def dshift(poly, j):
    return s.expand(sum(T[i + j] * s.diff(poly, T[i]) for i in range(1, 12)))


def source_l(poly):
    return s.expand((dshift(dshift(poly, 1), 1) - dshift(poly, 2)) / 2)


def source_h(poly):
    return s.expand(-sum(sum(s.binomial(i, l) * T[1 + l] * T[1 + i - l]
                            for l in range(1, i)) * s.diff(poly, T[i])
                        for i in range(1, 12)) / 2)


def source_a(n):
    return -sum(s.binomial(n, k) * T[1 + k] * T[1 + n - k] * T[n]
                for k in range(1, n))


R = {2: s.S.Zero}
for n in range(2, 11):
    R[n + 1] = s.expand(source_a(n) + source_l(R[n]) + source_h(R[n]))
    assert not any(R[n + 1].has(T[i]) for i in range(n + 1, 15))

# Fidelity: the full displayed R_2 through R_6, not only a guessed coefficient.
source_small = {
    2: 0,
    3: -2*T[2]**3,
    4: -12*T[2]*T[3]**2 + 6*T[2]**4,
    5: -20*T[2]*T[4]**2 - 30*T[3]**2*T[4] + 120*T[2]**2*T[3]**2 - 24*T[2]**5,
    6: -30*T[2]*T[5]**2 - 120*T[3]*T[4]*T[5] + 900*T[2]*T[3]**2*T[4]
       + 300*T[2]**2*T[4]**2 - 30*T[4]**3 - 1200*T[2]**3*T[3]**2
       + 210*T[3]**4 + 120*T[2]**6,
}
for n, stated in source_small.items():
    same(R[n], stated)
print('FIDELITY: full source R_2,...,R_6 reproduced from recurrence')
a = {n: s.expand(R[n + 1]).coeff(T[n], 1) for n in (6, 8, 10)}
same(a[8], -630*T[5]**2 - 504*T[3]*T[7] - 1008*T[4]*T[6])
assert len(s.Poly(R[9], *T).terms()) == 42
print('FIDELITY: Example 1.6 has 42 terms and its exact a_8 coefficient agrees')
# The source prints -210 T_3 T_6 in Remark 1.5; recurrence gives T_3 T_5.
same(a[6], -140*T[4]**2 - 210*T[3]*T[5])
print('SOURCE TYPO: Remark 1.5 prints T_3*T_6; recurrence gives T_3*T_5; both vanish on symmetric laws')

# Exact moments of masses p/2, 1-p, p/2 at -1,0,1 and source Prop. 2.6.
mu = {0: s.S.One}
for n in range(1, 11):
    mu[n] = s.expand((p/2)*(-1)**n + (1-p)*0**n + (p/2)*1**n)
kappa = {0: s.S.Zero}
for n in range(1, 11):
    kappa[n] = s.expand(mu[n] - sum(s.binomial(n-1, j-1)*kappa[j]*mu[n-j]
                                  for j in range(1, n)))
# Independent truncated formal log of M(t)=1-p+p*cosh(t).
u = sum(p*t**(2*j)/s.factorial(2*j) for j in range(1, 6))
log_m = s.Poly(s.expand(sum((-1)**(j+1)*u**j/Q(j) for j in range(1, 6))), t)
for n in range(1, 11):
    same(kappa[n], log_m.coeff_monomial(t**n)*s.factorial(n))
expected = {2: p, 4: p-3*p**2, 6: p-15*p**2+30*p**3,
            8: p-63*p**2+420*p**3-630*p**4}
for n, expr in expected.items():
    same(kappa[n], expr)
    print('K_%d = %s' % (n, s.factor(expr)))
assert all(kappa[n] == 0 for n in (1, 3, 5, 7, 9))
av = {n: s.expand(a[n].subs({T[j]: kappa[j] for j in range(2, 10)}, simultaneous=True))
      for n in a}
same(av[6], -140*kappa[4]**2)
same(av[8], -1008*kappa[4]*kappa[6])
same(av[10], -132*(21*kappa[6]**2+20*kappa[4]*kappa[8]))
P = 56700*p**4 - 56700*p**3 + 18165*p**2 - 1950*p + 41
same(21*kappa[6]**2 + 20*kappa[4]*kappa[8], p**2*P)
lo, hi = Q(1, 100), Q(1, 20)
assert P.subs(p, lo) == Q(23260367, 1000000) > 0
assert P.subs(p, hi) == -Q(28513, 1600) < 0
# Upper bound on P' throughout [lo,hi]: discard -170100*p^2;
# bound both positive monomials at hi.
same(s.diff(P, p), 226800*p**3-170100*p**2+36330*p-1950)
derivative_upper = 226800*hi**3 + 36330*hi - 1950
assert derivative_upper == -Q(2103, 20) < 0
assert s.Poly(P, p).count_roots(lo, hi) == 1
assert 0 < lo < hi < 1
assert 1-3*hi == Q(17, 20) > 0
assert 1-15*hi == Q(1, 4) > 0
for n in (2, 4, 6, 8):
    assert s.gcd(s.Poly(P, p), s.Poly(kappa[n], p)).degree() == 0
same(av[10] + 132*p**2*P)
print('P =', P)
print('P(1/100) =', P.subs(p, lo), '; P(1/20) =', P.subs(p, hi))
print("P' <=", derivative_upper, '< 0 on [1/100,1/20]; exact Sturm root count = 1')
print('At the unique root p0: K2,K4,K6 > 0; K8 = -21*K6^2/(20*K4) < 0')
print('At p0: a6 < 0, a8 < 0, a10 = 0; no zero-K4 or Gaussian loophole')
# Elementary moment determinacy of the three-point law: zero expectation of a
# nonnegative polynomial forces any law with these moments onto {-1,0,1}.
x = s.symbols('x', real=True)
nonnegative = x**2*(x**2-1)**2
same(nonnegative, x**6-2*x**4+x**2)
same(mu[6]-2*mu[4]+mu[2])
assert set(s.solve(nonnegative, x)) == {-1, 0, 1}
wminus, wzero, wplus = s.symbols('wminus wzero wplus')
assert s.solve([wminus+wzero+wplus-1, -wminus+wplus,
                wminus+wplus-p], (wminus,wzero,wplus)) == {
                    wminus:p/2, wzero:1-p, wplus:p/2}
print('Three-point moment-determinacy certificate: E[X^2*(X^2-1)^2]=0; masses uniquely recovered')

# Gaussian addition: its log-MGF adds sigma^2*t^2/2 exactly. This checks all
# orders needed for the counterexample without a numerical approximation.
smooth_log = log_m.as_expr()+sigma**2*t**2/2
for n in range(3, 11):
    same(s.expand(smooth_log).coeff(t,n)*s.factorial(n), kappa[n])
same(s.expand(smooth_log).coeff(t,2)*2, p+sigma**2)
for n, degree in [(6,8),(8,10),(10,12)]:
    scaled = a[n].subs({T[j]: b**j*kappa[j] for j in range(2,10)}, simultaneous=True)
    same(scaled, b**degree*av[n])
print('Gaussian addition preserves cumulants of orders 3..10 and the cancellation')
print('Scaling multiplies a6,a8,a10 by b^8,b^10,b^12 respectively')
# Symbolic Gaussian-mixture log-Hessian in posterior weights.
y, v = s.symbols('y v', real=True)
wm, w0, wp = s.symbols('wm w0 wp', real=True)
f0 = wm+w0+wp
f1 = (wm*(-1-y) + w0*(-y) + wp*(1-y))/v
f2 = wm*((-1-y)**2/v**2-1/v) + w0*(y**2/v**2-1/v) + wp*((1-y)**2/v**2-1/v)
posterior_var = (wm+wp)/f0 - ((wp-wm)/f0)**2
assert s.cancel(f2/f0-(f1/f0)**2 - (-1/v+posterior_var/v**2)) == 0
same((1-posterior_var).subs(w0,1-wm-wp), 1-wm-wp+(wp-wm)**2)
# With normalized nonnegative weights, 1-Var = w0+(wp-wm)^2 >= 0.
# Therefore sigma>1 gives (log f)'' <= -(sigma^2-1)/sigma^4 < 0.
print("Log-density Hessian = -1/sigma^2 + posterior_Var/sigma^4")
print('Posterior support [-1,1]: 1-Var = w0+(wp-wm)^2 >= 0')
print("sigma>1: (log f)'' <= -(sigma^2-1)/sigma^4 < 0")
normalized_k4 = kappa[4]/(p+sigma**2)**2
assert s.cancel(s.diff(normalized_k4,sigma) + 4*sigma*kappa[4]/(p+sigma**2)**3) == 0
print('Normalized K4 is strictly decreasing for sigma>1, proving distinct normalized laws')
print('Analytic obligations beyond exact algebra: finite Gaussian-mixture density is positive analytic;')
print('its everywhere finite MGF implies moment determinacy (classical MGF criterion).')
print('ALL_CHECKS_OK')
