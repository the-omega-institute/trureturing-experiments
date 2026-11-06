#!/usr/bin/env python3
"""Exact consistency checks, not a proof of an infinite-dimensional spectrum."""
from fractions import Fraction
import sympy as s

u = s.symbols('u', positive=True)
t, d, xi = s.symbols('t d xi', real=True)
R = s.symbols('R', positive=True)
half = s.Rational(1, 2)
etas = [s.Rational(n, 2) for n in range(1, 9)]

def same(a, b):
    delta = s.simplify(s.expand_func(a - b))
    assert delta == 0, (a, b, delta)

print('SymPy version:', s.__version__)
print('Arithmetic: exact integers, rationals, symbolic pi/Gamma/exponentials; no Float inputs')

# Fidelity: source equation (14.2), the first two conditional moments.
for eta in etas:
    mgf = (s.cos(t / eta) - d * s.sin(t / eta)) ** (-eta)
    mean = s.diff(mgf, t).subs(t, 0)
    second = s.diff(mgf, t, 2).subs(t, 0)
    same(mean, d)
    same(second - mean**2, (1 + d**2) / eta)
print('FIDELITY conditional mean and variance: eta = 1/2, 1, ..., 4 (8 exact cases)')

# Fidelity: limiting Cauchy interval-mass ratio in Proposition 14.7.
ratio = (s.atan(3 * R / 2) - s.atan(R)) / (s.atan(2 * R) - s.atan(R))
same(s.limit(ratio, R, s.oo), s.Rational(2, 3))
print('FIDELITY Proposition 14.7 interval mass ratio: 2/3')

# The gamma edge density is normalized and has mean one.
for eta in etas:
    density = eta**eta / s.gamma(eta) * u**(eta - 1) * s.exp(-eta*u)
    moments = [s.integrate(u**j * density, (u, 0, s.oo)) for j in range(3)]
    same(moments[0], 1)
    same(moments[1], 1)
    same(moments[2] - 1, 1/eta)
print('FIDELITY edge gamma normalization, mean, variance: 8 exact cases')

# Expected settlement: after u = exp(v), the actual Mellin/Fourier integral
# has exponent eta - 1/2 + i xi. xi remains an arbitrary real symbol.
for eta in etas:
    integrand = eta**eta / s.gamma(eta) * u**(eta - half + s.I*xi) * s.exp(-eta*u)
    actual = s.integrate(integrand, (u, 0, s.oo))
    expected = s.gamma(eta + half + s.I*xi) / (s.gamma(eta)*eta**(half+s.I*xi))
    assert not actual.has(s.Integral, s.Float), actual
    same(actual, expected)
    print('SYMBOL eta =', eta, ': Mellin integral equals proposed Gamma symbol for symbolic real xi')

# The same identity at all small nonnegative Mellin moments.
for eta in etas:
    for j in range(7):
        actual = eta**eta / s.gamma(eta) * s.integrate(
            u**(eta-half+j)*s.exp(-eta*u), (u, 0, s.oo))
        expected = s.gamma(eta+half+j)/(s.gamma(eta)*eta**(half+j))
        same(actual, expected)
print('SYMBOL nonnegative Mellin moments: 8 x 7 = 56 exact cases')

# Jacobian/normalization audit of the logarithmic limiting kernel.
v = s.symbols('v', real=True)
eta = s.symbols('eta', positive=True)
x = s.symbols('x', positive=True)
k = eta**eta / s.gamma(eta) * s.exp((eta+half)*v-eta*s.exp(v))
edge = eta**eta/s.gamma(eta)*u**(eta-1)*s.exp(-eta*u)
same(k.subs(v, s.log(u))/u, s.sqrt(u)*edge)
same(s.integrate(u**(2*eta)*s.exp(-2*eta*u), (u, 0, s.oo)),
     s.gamma(2*eta+1)/(2*eta)**(2*eta+1))
print('LOG KERNEL: Jacobian gives sqrt(u) times edge density; squared-kernel integral finite')

# Analytic source of the s-tail O(exp(-2s)) error, using z = exp(-s).
z = s.symbols('z', positive=True)
same(s.limit((1-s.atan(z)/z)/z**2, z, 0), s.Rational(1, 3))
same(s.limit(((1+z**2)**(-(eta+1)/2)-1)/z**2, z, 0), -(eta+1)/2)
print('REMAINDER elementary coefficients: a(s) = 1 - exp(-2s)/3 + ..., beta(s) = 1 - (eta+1)exp(-2s)/2 + ...')

# Fredholm parametrix scalar identity. This checks the sign and the
# necessary constant term; an inverse symbol is not itself Schwartz.
w, lam = s.symbols('w lam', nonzero=True)
same(1/(w-lam), -1/lam + w/(lam*(w-lam)))
print('PARAMETRIX: inverse symbol = -1/lambda + phi/[lambda(phi-lambda)]')

# Exact checks of the proposed norm constants and strict bound in the
# three nominated cases, using only the rigorous rational bounds 3 < pi < 22/7.
expected_squares = {half: 2/s.pi, s.Integer(1): s.pi/4, s.Integer(3): 75*s.pi/256}
for a, expected in expected_squares.items():
    M = s.gamma(a+half)/(s.sqrt(a)*s.gamma(a))
    same(M**2, expected)
assert Fraction(2, 3) < 1
assert Fraction(22, 7) / 4 < 1
assert Fraction(75, 256) * Fraction(22, 7) < 1
print('NORM eta = 1/2, 1, 3: exact squares 2/pi, pi/4, 75pi/256 are strictly below one')
# Quantitative Gamma remainder normalization. Brent's log-Gamma theorem
# supplies the analytic O(|z|^-3) error; here only the exact algebra that
# cancels an apparent 1/x term and yields the 1/x^2 coefficient is checked.
a = s.symbols('a', positive=True)
zeta = s.symbols('zeta', positive=True)
log_modulus_correction = ((2*a-1)/2*s.log(1+a*a*zeta*zeta)
    + 2/zeta*s.atan(a*zeta)-2*a + a*zeta*zeta/(6*(1+a*a*zeta*zeta)))
coefficient_b = s.limit(log_modulus_correction/zeta**2, zeta, 0)
coefficient_x = s.simplify(4*coefficient_b.subs(a, eta/2))
same(coefficient_x, eta*(eta-1)*(eta-2)/6)
same(s.Rational(1,2)+coefficient_x/eta**2, (eta**2+2)/(6*eta))
print('GAMMA REMAINDER algebra: density coefficient eta(eta-1)(eta-2)/6; alpha coefficient (eta^2+2)/(6eta)')

# Exact numerical-free bounds appearing in the convention-separation
# witness eta=1, lambda=sqrt(pi)/4. The Fredholm index theorem itself is
# not certified by these scalar inequalities.
assert sum(Fraction(6**(2*j), s.factorial(2*j)) for j in range(3)) == 73
assert 73 > 68
assert Fraction(220,81) < 3
print('CONVENTION GUARD scalar bounds: cosh(6)>73>68 and 220/81<3<pi')

print('LIMIT: exact checks test fidelity and symbol identities, not compactness or Fredholmness')
print('ALL_CHECKS_OK')
