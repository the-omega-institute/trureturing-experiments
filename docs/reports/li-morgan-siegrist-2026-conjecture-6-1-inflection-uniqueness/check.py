"""Li-Morgan-Siegrist, J. Anal. (2026), Conjecture 6.1: exact counterexamples.

For alpha = 3/2 and a polynomial forcing f, the solution of (7) is
F(x) = I^{3/2} f(x) - I^{3/2} f(1) sqrt(x).  For f = sum_j c_j x^j,
I^{3/2} x^j = Gamma(3/2)^{-1} int_0^x (x - t)^{1/2} t^j dt is computed here by the
binomial expansion of t^j = (x - u)^j and int_0^x u^{i + 1/2} du (no Gamma closed form),
and checked against the Beta formula Gamma(j + 1) / Gamma(j + 5/2) x^{j + 3/2}.
Then P(x) = Gamma(3/2) x^{3/2} F''(x) is an exact rational polynomial, and its signs
at rational points certify sign-changing zeros of F'' in (0, 1).
"""
import sympy as sp
import mpmath as mp

x, u, A = sp.symbols('x u A', positive=True)
HALF = sp.Rational(1, 2)
G32 = sp.gamma(sp.Rational(3, 2))


def frac_int_monomial(j):
    # Gamma(3/2) * I^{3/2} x^j = int_0^x u^{1/2} (x - u)^j du
    s = sum(sp.binomial(j, i) * (-1) ** i * x ** (j - i) * x ** (i + sp.Rational(3, 2)) / (i + sp.Rational(3, 2))
            for i in range(j + 1))
    return sp.simplify(s)


def curvature_poly(f):
    cs = sp.Poly(f, x).all_coeffs()[::-1]
    g = sum(c * frac_int_monomial(j) for j, c in enumerate(cs))  # Gamma(3/2) * I^{3/2} f
    for j in range(len(cs)):
        beta = sp.gamma(j + 1) / sp.gamma(j + sp.Rational(5, 2)) * x ** (j + sp.Rational(3, 2)) * G32
        assert sp.simplify(frac_int_monomial(j) - beta) == 0, j
    gF = g - g.subs(x, 1) * sp.sqrt(x)  # Gamma(3/2) * F
    P = sp.expand(sp.simplify(x ** sp.Rational(3, 2) * sp.diff(gF, x, 2)))
    return sp.Poly(P, x)


def hyps(f):
    # polynomial, so C([0,1]) and C^1((0,1]) hold, and x^3 f'(x) -> 0; check f(1)=0, f>=0, f not 0
    assert f.subs(x, 1) == 0
    q = sp.Poly(sp.cancel(f / (1 - x)), x)
    assert all(c >= 0 for c in q.all_coeffs()) and q.eval(0) > 0  # f = (1-x) q(x) with q > 0 on [0,1]
    assert sp.limit(x ** 3 * sp.diff(f, x), x, 0, '+') == 0


ok = True
# Instance used for the Lean refutation
f = sp.expand((1 - x) * (1 + 5 * x ** 12)); hyps(f)
P = curvature_poly(f)
print('instance f = (1 - x)(1 + 5 x^12); P(x) =', P.as_expr())
pts = [sp.Rational(7, 10), sp.Rational(3, 4), sp.Rational(4, 5), sp.Rational(9, 10)]
sg = [sp.sign(P.eval(p)) for p in pts]
print('signs at 7/10, 3/4, 4/5, 9/10:', sg, [str(P.eval(p)) for p in pts])
ok &= sg == [1, -1, 1, -1]
rr = [r for r in P.real_roots() if 0 < r < 1]
print('zeros of F\'\' in (0,1):', [sp.N(r, 12) for r in rr], 'all simple:', all(P.diff().eval(r) != 0 for r in rr))
ok &= len(rr) == 3
# Independent quadrature check of F'' (no polynomial formula) at three points
mp.mp.dps = 40
fn = sp.lambdify(x, f, 'mpmath')
I = lambda xx: mp.quad(lambda t: (xx - t) ** mp.mpf('0.5') * fn(t), [0, xx]) / mp.gamma(1.5)
C1 = I(mp.mpf(1)); F = lambda xx: I(xx) - C1 * mp.sqrt(xx)
for xx in ['0.7', '0.75', '0.8', '0.9']:
    d2 = mp.diff(F, mp.mpf(xx), 2) * mp.gamma(1.5) * mp.mpf(xx) ** 1.5
    err = abs(d2 - P.eval(sp.Rational(xx)))
    print('quadrature x =', xx, 'P =', mp.nstr(d2, 15), 'abs err', mp.nstr(err, 3)); ok &= err < mp.mpf('1e-15')
# Search seat's family: f_A = (1 - x)(1 + A x^20), A in [19, 21]; P is affine in A
fA = sp.expand((1 - x) * (1 + A * x ** 20))
PA = curvature_poly(fA)
p = sp.Rational(2 ** 39) / sp.binomial(40, 20)
claimed = HALF / 5 + x / 2 - x ** 2 + A * p * (sp.Rational(1, 26445) + x ** 21 * (1 - sp.Rational(42, 41) * x))
print('family P_A matches the closed form:', sp.expand(PA.as_expr() - claimed) == 0); ok &= sp.expand(PA.as_expr() - claimed) == 0
for a in (19, 21):
    hyps(fA.subs(A, a))
    s = [sp.sign(PA.as_expr().subs({A: a, x: v})) for v in (HALF, sp.Rational(3, 4), sp.Rational(19, 20), 1)]
    print('A =', a, 'signs at 1/2, 3/4, 19/20, 1:', s); ok &= s == [1, -1, 1, -1]
# Other orders: same shape (1 - x)(1 + A x^k); x^{3-alpha} F'' is a polynomial with Gamma-ratio coefficients
def sign_changes(al, k, a, n=4000):
    al = mp.mpf(al); c = {0: 1, 1: -1, k: a, k + 1: -a}
    Cn = sum(v * mp.gamma(j + 1) / mp.gamma(j + al + 1) for j, v in c.items())
    Q = lambda t: sum(v * mp.gamma(j + 1) / mp.gamma(j + al - 1) * t ** (j + 1) for j, v in c.items()) - Cn * (al - 1) * (al - 2)
    vals = [Q(mp.mpf(i) / n) for i in range(1, n)]
    return sum(1 for a1, b1 in zip(vals, vals[1:]) if a1 * b1 < 0)
for al, k, a in [('1.1', 8, 10), ('1.3', 8, 5), ('1.7', 30, 10)]:
    n = sign_changes(al, k, a); print('alpha =', al, 'k =', k, 'A =', a, 'sign changes on a 1/4000 grid:', n); ok &= n == 3
print('RESULT', 'PASS' if ok else 'FAIL')
