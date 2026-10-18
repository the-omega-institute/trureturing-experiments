"""Maity-Furuichi-Dutta, arXiv:2607.01805v1, Eq. (9) and Section 4, Case 5: for non-commuting density operators,
D_alpha(rho||sigma) <= B1(rho||sigma) <= B2(rho||sigma), with B1 = D^(r)_{alpha,z} (Definition 2) and
B2 = -Tr[rho ln_{1-alpha}((rho^{-q/2} sigma^q rho^{-q/2})^{1/q})] (Eq. (47)), for 0 <= alpha < 1 and q > 1 - alpha.
Counterexample (exact): alpha = r = 1/2, z = 2, q = 1, rho = diag(4/5, 1/5), sigma = (1/125)[[28, 6 sqrt 6], [6 sqrt 6, 97]].
The program evaluates B1, B2 and D_alpha with sympy exactly (closed forms) and with mpmath at 50 digits, and checks
B1 > 2/5 > B2 > D_alpha = 48/125. It also reports how closely the same code reproduces the paper's four Case 3 numbers (reported, not a pass condition).
Requires sympy and mpmath.
"""
import sympy as sp, mpmath as mp
mp.mp.dps = 50
s6 = sp.sqrt(6)
rho = sp.Matrix([[sp.Rational(4,5),0],[0,sp.Rational(1,5)]])
sig = sp.Matrix([[28,6*s6],[6*s6,97]])/125
def mpow(M, t):
    # spectral power of a symmetric positive definite matrix, high precision
    A = mp.matrix([[mp.mpf(sp.N(M[i,j],60)) for j in range(M.shape[1])] for i in range(M.shape[0])])
    E, Q = mp.eighe(A)
    return Q*mp.diag([e**t for e in E])*Q.T
def mpowm(A, t):
    E, Q = mp.eighe(A)
    return Q*mp.diag([e**t for e in E])*Q.T
def lnr(x, r): return (x**r - 1)/r
def B1(rho, sig, a, r, z):
    s = mpow(sig, (1-a)/(2*z)); m = s*mpow(rho, a/z)*s
    F = sum(mpowm(m, z)[i,i] for i in range(2))
    return lnr(F, r)/(a-1)
def B2(rho, sig, a, q):
    R = mpow(rho, -q/2); A = R*mpow(sig, q)*R
    L = mpowm(A, 1/q)
    E, Q = mp.eighe(L); Lg = Q*mp.diag([lnr(e, 1-a) for e in E])*Q.T
    P = mpow(rho, 1)
    return -sum((P*Lg)[i,i] for i in range(2))
def Dts(rho, sig, a):
    X = mpow(rho, a)*mpow(sig, 1-a)
    return (1 - sum(X[i,i] for i in range(2)))/(1-a)
a, r, z, q = mp.mpf(1)/2, mp.mpf(1)/2, 2, 1
b1, b2, d = B1(rho, sig, a, r, z), B2(rho, sig, a, q), Dts(rho, sig, a)
print("commutator nonzero:", (rho*sig - sig*rho) != sp.zeros(2))
print("tr rho, tr sigma, det sigma:", rho.trace(), sp.nsimplify(sig.trace()), sp.simplify(sig.det()))
print("B1 =", mp.nstr(b1, 20), " closed form 4(1 - sqrt((2117+288 sqrt2)/3125)) =", mp.nstr(4*(1-mp.sqrt((2117+288*mp.sqrt(2))/3125)), 20))
print("B2 =", mp.nstr(b2, 20), " closed form 2 - 20/sqrt(154) =", mp.nstr(2-20/mp.sqrt(154), 20))
print("D_alpha =", mp.nstr(d, 20), " (48/125 =", mp.nstr(mp.mpf(48)/125, 20), ")")
ok = b1 > mp.mpf(2)/5 > b2 > d and abs(d - mp.mpf(48)/125) < mp.mpf(10)**-40
# reproduce the paper's Case 3 numbers (rho = [[0.6,0.3],[0.3,0.4]], sigma = [[0.4,0.2],[0.2,0.4]])
r3 = sp.Matrix([[sp.Rational(6,10),sp.Rational(3,10)],[sp.Rational(3,10),sp.Rational(4,10)]])
s3 = sp.Matrix([[sp.Rational(4,10),sp.Rational(2,10)],[sp.Rational(2,10),sp.Rational(4,10)]])
c3 = [(0.01, 0.207025, 0.203554), (0.99, 0.252861, 0.264951)]
for al, p1, p2 in c3:
    x1 = B1(r3, s3, mp.mpf(al), mp.mpf('0.7'), 2); x2 = B2(r3, s3, mp.mpf(al), mp.mpf('1.5'))
    print("Case 3 alpha=%s: B1=%s (paper %s) B2=%s (paper %s)" % (al, mp.nstr(x1, 7), p1, mp.nstr(x2, 7), p2))
    print('  agreement within 5e-6: B1', abs(x1 - p1) < 5e-6, 'B2', abs(x2 - p2) < 5e-6)
# Mechanism: Theorem 2 gives B2 < B1 for commuting rho != sigma; a small rotation keeps the strict order.
def rot(th, d1, d2):
    c, s = sp.cos(th), sp.sin(th)
    R = sp.Matrix([[c, -s], [s, c]]); return R*sp.diag(d1, d2)*R.T
near = 0
grid = [(sp.Rational(7,10), sp.Rational(4,10)), (sp.Rational(9,10), sp.Rational(1,2)), (sp.Rational(6,10), sp.Rational(2,10))]
for d1, e1 in grid:
    for th in (sp.Rational(1,100), sp.Rational(1,10)):
        r_ = sp.diag(d1, 1-d1); s_ = rot(th, e1, 1-e1)
        x1 = B1(r_, s_, mp.mpf(1)/2, mp.mpf(1)/2, 2); x2 = B2(r_, s_, mp.mpf(1)/2, 1)
        near += x1 > x2
print("near-commuting pairs with B1 > B2: %d of %d" % (near, 2*len(grid)))
ok = ok and near == 2*len(grid)
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
