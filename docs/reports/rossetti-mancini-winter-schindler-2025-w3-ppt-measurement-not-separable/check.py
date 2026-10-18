#!/usr/bin/env python3
"""Exact check: the optimal PPT effect pair (M', 1 - M') of Rossetti-Mancini-Winter-Schindler
(arXiv:2510.10058v1, Sec. X.0.4, Eq. (83)) is not fully separable, although both effects are PPT.
M' = Q1 + Q2/2 + (2/3) Q5 + (1/12) Q6; R = 1 - M'. Exit 0 iff every check passes."""
import itertools, sys
import sympy as sp

def ket(bits):
    v = sp.zeros(8, 1); v[int(bits, 2)] = 1; return v
def proj(v):
    return v * v.H
I8 = sp.eye(8)
W = (ket('100') + ket('010') + ket('001')) / sp.sqrt(3)
Wb = (ket('011') + ket('101') + ket('110')) / sp.sqrt(3)
Q1, Q2, Q5, Q6 = proj(W), proj(Wb), proj(ket('000')), proj(ket('111'))
Mp = sp.simplify(Q1 + Q2 / 2 + sp.Rational(2, 3) * Q5 + sp.Rational(1, 12) * Q6)
R = sp.simplify(I8 - Mp)
fails = []
def check(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)
    if not ok: fails.append(name)

check('trace M\' = 9/4', sp.simplify(Mp.trace() - sp.Rational(9, 4)) == 0)
Rshown = sp.Matrix([[4,0,0,0,0,0,0,0],[0,8,-4,0,-4,0,0,0],[0,-4,8,0,-4,0,0,0],[0,0,0,10,0,-2,-2,0],
                    [0,-4,-4,0,8,0,0,0],[0,0,0,-2,0,10,-2,0],[0,0,0,-2,0,-2,10,0],[0,0,0,0,0,0,0,11]]) / 12
check('R equals the displayed matrix', sp.simplify(R - Rshown) == sp.zeros(8, 8))

def ptranspose(X, i):
    Y = sp.zeros(8, 8)
    for r, c in itertools.product(range(8), repeat=2):
        rb, cb = list(format(r, '03b')), list(format(c, '03b'))
        rb[i], cb[i] = cb[i], rb[i]
        Y[int(''.join(rb), 2), int(''.join(cb), 2)] = X[r, c]
    return Y
def psd(X):
    ev = X.eigenvals()
    return all(sp.re(sp.nsimplify(e)) >= 0 and sp.simplify(sp.im(e)) == 0 for e in ev), ev
for name, X in (('M\'', Mp), ('R', R)):
    ok, ev = psd(X); check(f'{name} >= 0 (eigenvalues {sorted(ev, key=lambda e: float(e))})', ok)
    for i in range(3):
        ok, ev = psd(ptranspose(X, i)); check(f'{name}^T{i+1} >= 0 (eigenvalues {sorted(ev, key=lambda e: float(e))})', ok)

w = ket('100') + ket('010') + ket('001')
z = [2*ket('000') + ket('110') + ket('101'), 2*ket('000') + ket('110') + ket('011'), 2*ket('000') + ket('101') + ket('011')]
check('R w = 0', sp.simplify(R * w) == sp.zeros(8, 1))
for i in range(3):
    check(f'R^T{i+1} z{i+1} = 0', sp.simplify(ptranspose(R, i) * z[i]) == sp.zeros(8, 1))
P1 = proj(ket('100')) + proj(ket('010')) + proj(ket('001'))
P2 = proj(ket('011')) + proj(ket('101')) + proj(ket('110'))
t1, t2 = sp.simplify((P1 * R).trace()), sp.simplify((P2 * R).trace())
check(f'Tr(P1 R) = 2 (got {t1})', t1 == 2)
check(f'Tr(P2 R) = 5/2 (got {t2})', t2 == sp.Rational(5, 2))
check('Tr(P2 R) != 2 Tr(P1 R)', t2 != 2 * t1)

# Product-vector lemma: the four kernel conditions force ||P2 u||^2 = 2 ||P1 u||^2.
a0, a1, b0, b1, c0, c1 = sp.symbols('a0 a1 b0 b1 c0 c1')
def kron3(a, b, c):
    return sp.Matrix([a[i] * b[j] * c[k] for i in range(2) for j in range(2) for k in range(2)])
u = kron3([a0, a1], [b0, b1], [c0, c1])
conjs = [kron3([sp.conjugate(a0), sp.conjugate(a1)], [b0, b1], [c0, c1]),
         kron3([a0, a1], [sp.conjugate(b0), sp.conjugate(b1)], [c0, c1]),
         kron3([a0, a1], [b0, b1], [sp.conjugate(c0), sp.conjugate(c1)])]
print('INFO <w|u> =', sp.expand((w.T * u)[0]))
for i in range(3):
    print(f'INFO <z{i+1}|u^*{i+1}> =', sp.expand((z[i].T * conjs[i])[0]))
# Concrete solutions with a0 b0 c0 = 1: x + y + z = 0 and |x|^2 = |y|^2 = |z|^2 = 2.
om = sp.Rational(-1, 2) + sp.sqrt(3) * sp.I / 2
for (x, y, zz) in ((sp.sqrt(2), sp.sqrt(2) * om, sp.sqrt(2) * om**2), (sp.sqrt(2) * sp.I, sp.sqrt(2) * sp.I * om, sp.sqrt(2) * sp.I * om**2)):
    sub = {a0: 1, b0: 1, c0: 1, a1: x, b1: y, c1: zz}
    vals = [sp.expand(sp.expand_complex((w.T * u)[0].subs(sub)))] + [sp.expand(sp.expand_complex((z[i].T * conjs[i])[0].subs(sub))) for i in range(3)]
    uu = u.subs(sub); n1 = sp.expand(sp.expand_complex((uu.H * P1 * uu)[0])); n2 = sp.expand(sp.expand_complex((uu.H * P2 * uu)[0]))
    check(f'sample solution satisfies the four kernel equations and ||P2u||^2 = 2||P1u||^2 ({n1}, {n2})', all(v == 0 for v in vals) and sp.simplify(n2 - 2 * n1) == 0)

# M' is fully separable: M' = (1/48) sum_k |u_k><u_k|, u_k = (sqrt2|0> + i^k|1>)^{(x)3}.
S = sp.zeros(8, 8)
for k in range(4):
    v = sp.Matrix([sp.sqrt(2), sp.I**k]); uk = kron3(v, v, v); S += uk * uk.H
check('M\' = (1/48) sum_k |u_k><u_k| (fully separable)', sp.simplify(S / 48 - Mp) == sp.zeros(8, 8))
print('RESULT', 'FAIL' if fails else 'PASS', len(fails))
sys.exit(1 if fails else 0)
