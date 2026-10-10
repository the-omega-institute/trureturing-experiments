# Check of the refutation of Fucci–Stanfill (arXiv:2411.17860v1, Remark B.3; Ann. Henri Poincaré 27 (2026) 2073–2116):
# "we in fact conjecture that all of the g_{m,j}(p/q,beta) are nonzero".  Claim: g_{q,1}(p/q,beta) = 0 for coprime 0<p<q.
import sympy as sp
from math import gcd
y, T, a_, x = sp.symbols('y T a x')
# general definitions (B.17, C_m, G_k via generalized Bernoulli polynomials, P_k, script-E_j)
def C(m, y):
    return sum(sp.binomial(2*m-1, 2*j) * 2**(2*m) * y**(2*j) / (m-j) * sp.bernoulli(2*(m-j)) for j in range(m))
def E(k, y):
    return 2*(2*y)**(2*k-1) - (2*y)**(2*k)/k - C(k, y)
t = sp.symbols('t')
def genB(n, a, xx):
    # (t/(e^t-1))^a e^{x t} = sum B_n^{(a)}(x) t^n/n!   (DLMF 24.16.1), via series
    F = sp.series(t/(sp.exp(t)-1), t, 0, n+2).removeO()
    u = sp.expand(F - 1)
    powa = sum(sp.binomial(a, k) * u**k for k in range(n+2))
    ser = sp.series(sp.expand(powa) * sp.exp(xx*t), t, 0, n+1).removeO()
    return sp.factorial(n) * sp.expand(ser).coeff(t, n)
def G(k, xx, yy):
    return sp.binomial(xx-yy, k) * genB(k, xx-yy+1, xx)   # G_k(x,y) = binom(x-y,k) B_k^{(x-y+1)}(x)
def P(k, yy):
    if k == 0: return sp.Integer(1)
    tot = 0
    for j in range(k+1):
        scrE = 1 if j == 0 else T/2 * E(j, yy)
        tot += 4**(k-j) * scrE * G(2*(k-j), 1-yy, yy)
    return sp.expand(tot)
P1 = sp.simplify(sp.expand_func(P(1, y)))
A_print = sp.Rational(2,3)*y*(2*y**2-3*y+1); B_print = (-12*y**2+12*y+2)/sp.Integer(6)
A_def = sp.simplify(P1.coeff(T, 0)); B_def = sp.simplify(P1.coeff(T, 1))
print('P1 constant part from the definitions equals B.21:', sp.simplify(A_def - A_print) == 0)
print('P1 T-coefficient from the definitions:', sp.factor(B_def), '; B.21 prints', sp.factor(B_print), '; difference', sp.simplify(B_def - B_print))
# both forms are symmetric under y -> 1-y
for Bf in (B_def, B_print):
    assert sp.simplify(Bf.subs(y, 1-y) - Bf) == 0
print('T-coefficient symmetric under y -> 1 - y (both forms)')
# g_{q,1} via the formal logarithm, with arbitrary higher Omega_k
results = []
for q in range(2, 9):
    for p in range(1, q):
        if gcd(p, q) != 1: continue
        nu = sp.Rational(p, q); yp = (1+nu)/2; ym = (1-nu)/2
        Cc = sp.Symbol('Cc')    # Omega_0 = 2^{2nu-1} Gamma(1+nu)/Gamma(-nu) cot(beta) e^{i pi nu}, nonzero for beta != pi/2
        Om = {0: Cc, 1: sp.expand(Cc*(P1.subs(y, yp)) - P1.subs(y, ym)*Cc)}
        K = 4
        for k in range(2, K+1):
            Om[k] = sum(sp.Symbol(f'w{k}_{i}') * T**i for i in range(k+1))   # arbitrary higher coefficients
        order = q + p
        def trunc(e):
            e = sp.expand(e)
            return sum(e.coeff(x, d) * x**d for d in range(order+1))
        W = trunc(sum(Om[k] * x**(p + k*q) for k in range(K+1) if p + k*q <= order))
        logser = 0; Wn = 1
        for n in range(1, order//p + 1):
            Wn = trunc(Wn * W)
            if Wn == 0: break
            logser += sp.Rational((-1)**(n+1), n) * Wn
        Sq = sp.expand(logser).coeff(x, order)
        g1 = sp.simplify(Sq.coeff(T, 1)); g0 = sp.simplify(Sq.coeff(T, 0))
        assert g1 == 0, (p, q, g1)
        assert sp.simplify(g0 - (Cc*nu*(nu**2-1)/3 if p > 1 else g0)) == 0
        assert not any(s.name.startswith('w') for s in Sq.free_symbols)
        results.append((p, q, sp.factor(g0)))
print('g_{q,1} = 0 for all coprime 0<p<q<=8 (higher Omega_k arbitrary); g_{q,0} values:', results[:6], '...')
print('ALL_OK')
