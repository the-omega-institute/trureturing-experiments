#!/usr/bin/env python3
"""Faust-Friedman-O'Malley-Ramos-Sharma (arXiv:2602.08304v1), Section 5, Conjecture 3, odd clause: "When q is odd, then
V(I') has only regular nonzero points".

Exact (sympy) check at q = 15 (m = 7). L_V is the Floquet matrix at z = 1 (diagonal V, cyclic off-diagonal entries 1),
D_V(lambda) = det(L_V - lambda I), zeta_U = D_U - D_0 for U = (u_1..u_7, 0, -u_7..-u_1), p'_k = (-1)^(n+k)[lambda^(n-2k)]
zeta_U. The point u = (a, b, 0, -b, -a, a, b) with a^2 = 1 - 2i, b = i a (the period-5 word (a, b, 0, -b, -a) repeated):
(1) the transfer-matrix identity det(lambda I - L_V) = tr(T(V_n)...T(V_1)) - 2 agrees with a direct determinant;
(2) every odd-offset coefficient of zeta_U vanishes (Section 4);
(3) p'_k(u) = 0 for k = 1..7 exactly;
(4) columns u_1 and u_6 of the Jacobian (d p'_k / d u_j)(u) are equal exactly, so its determinant is 0 (numerical rank 2).
Exit 0 and RESULT PASS when every check holds.
"""
import sympy as sp, random, sys
random.seed(2602)
lam = sp.Symbol('lam'); m = 7; n = 2 * m + 1
u = sp.symbols('u1:%d' % (m + 1)); V = list(u) + [0] + [-x for x in reversed(u)]
def T(v): return sp.Matrix([[lam - v, -1], [1, 0]])
def dl(Vs):
    M = sp.eye(2)
    for v in Vs: M = T(v) * M
    return sp.expand(M.trace() - 2)
P = dl(V); P0 = dl([0] * n)
vals = {x: sp.Rational(random.randint(-5, 5), random.randint(1, 4)) for x in u}; lv = sp.Rational(3, 7)
Ln = sp.zeros(n, n)
for i in range(n):
    Ln[i, i] = sp.sympify(V[i]).subs(vals); Ln[i, (i + 1) % n] += 1; Ln[(i + 1) % n, i] += 1
ok1 = sp.simplify((lv * sp.eye(n) - Ln).det() - P.subs(vals).subs(lam, lv)) == 0
zeta = sp.Poly(sp.expand((-1) ** n * (P - P0)), lam)
ok2 = all(sp.expand(zeta.coeff_monomial(lam ** (n - 2 * k - 1))) == 0 for k in range(0, m))
pk = [sp.expand((-1) ** (n + k) * zeta.coeff_monomial(lam ** (n - 2 * k))) for k in range(1, m + 1)]
a = sp.sqrt(sp.sqrt(5) / 2 + sp.Rational(1, 2)) - sp.I * sp.sqrt(sp.sqrt(5) / 2 - sp.Rational(1, 2)); b = sp.I * a
ok_a = sp.simplify(sp.expand(a ** 2) - (1 - 2 * sp.I)) == 0
pt = dict(zip(u, [a, b, 0, -b, -a, a, b]))
ok3 = all(sp.simplify(sp.expand(p.subs(pt))) == 0 for p in pk)
J = sp.Matrix([[sp.diff(p, x) for x in u] for p in pk]).subs(pt)
ok4 = all(sp.simplify(sp.expand(J[i, 0] - J[i, 5])) == 0 for i in range(m))
rank = sp.Matrix(sp.N(J, 50)).rank(iszerofunc=lambda x: abs(x) < 1e-30)
print('(1) transfer identity:', ok1); print('(2) odd coefficients vanish:', ok2); print('    a^2 = 1 - 2i:', ok_a)
print('(3) p_k(u) = 0 for k = 1..7:', ok3); print('(4) Jacobian columns 1 and 6 equal:', ok4, '; numerical rank', rank)
ok = ok1 and ok2 and ok_a and ok3 and ok4 and rank < m
print('RESULT PASS' if ok else 'RESULT FAIL'); sys.exit(0 if ok else 1)
