#!/usr/bin/env python3
"""Gibilisco–Isola math-ph/0407007v2, Conjecture 4.2, second clause at p = 3 (noncommutative alpha-geometry).

The alpha-geometry on D_n^1 (positive definite n x n density matrices) is the pull-back of the real
Hilbert–Schmidt inner product Re Tr(AB) under A_3(rho) = 3 rho^(1/3). With y_i = lambda_i^(1/3) and s_k = sum y_i^k
(so s_3 = 1), the proposed curvature formula is
    Scal_3 = N / (9 s_4^2),  N = [(4n^2 - 2) s_1^2 - 2n s_2] s_4 - 8n s_1 s_5 + 8 s_6.

(1) Intrinsic scalar curvature of the pull-back metric, computed in coordinates (rho = I/n + sum x_k G_k / 2 with a
    trace-orthogonal Hermitian basis G_k), from Christoffel symbols and the Ricci contraction by nested central
    differences at 60 digits; the embedding derivative uses the Daleckii–Krein formula. Points for n = 2 and n = 3,
    compared with N/(9 s_4^2). Control: the same code for the map 2 rho^(1/2) gives (n^2-1)(n^2-2)/4 (source).
(2) The Gauss-equation evaluation (tr S)^2 - tr(S^2) for the hypersurface Tr Y^3 = 1 in the real Hermitian space,
    with the shape operator built as an explicit n^2 x n^2 matrix, at random points for n = 2..5, compared with N/s_4^2.
(3) Symbolic algebra: the chain-rule identity d/dlambda_j R - d/dlambda_i R = 4(x-y) Q / (27 x^2 y^2 s_4^3) for
    n = 2, 3, 4; divisibility defining Q; the moment-deficit identity
    Q = L + 2n x^2 y^2 (c + x^4 + y^4) u + 2[8a(n-2) x^2 y^2 + n c (x+y) + n H_5(1,1,8)] v + 16 x^2 y^2 w
    with symbolic n; the displayed closed form of L; positivity of L for n >= 4 (n = k + 4), n = 2 and n = 3
    (after subtracting y^2 z^4 (2x^2 - z^2)^2).
(4) Strict increase of N/(9 s_4^2) under random nontrivial T-transforms (n = 2..7)."""
import sys, random, itertools
import mpmath as mp
import sympy as sp

bad = 0
def report(label, ok):
    global bad
    print(label, ok)
    bad += 0 if ok else 1

def N_of(ys, n):
    s = lambda k: sum(t**k for t in ys)
    return ((4*n*n - 2)*s(1)**2 - 2*n*s(2))*s(4) - 8*n*s(1)*s(5) + 8*s(6)

def pred_p3(lams):
    ys = [mp.cbrt(l) for l in lams]; n = len(ys)
    return N_of(ys, n) / (9 * sum(t**4 for t in ys)**2)

# ---------------- (1) intrinsic curvature in coordinates ----------------
mp.mp.dps = 60
def herm_basis(n):
    """Hermitian G_k, trace-orthogonal with Tr(G_k G_l) = 2 delta_kl, traceless (n^2 - 1 of them)."""
    B = []
    for i in range(n):
        for j in range(i + 1, n):
            M = mp.zeros(n); M[i, j] = 1; M[j, i] = 1; B.append(M)
            M = mp.zeros(n); M[i, j] = -1j; M[j, i] = 1j; B.append(M)
    for k in range(1, n):
        M = mp.zeros(n)
        c = mp.sqrt(mp.mpf(2) / (k * (k + 1)))
        for i in range(k): M[i, i] = c
        M[k, k] = -k * c
        B.append(M)
    return B

def metric(n, G, x, f, fp):
    rho = mp.eye(n) / n
    for k, g in enumerate(G): rho += x[k] * g / 2
    E, U = mp.eighe(rho)
    lam = [mp.re(E[i]) for i in range(n)]
    F = mp.matrix(n, n)
    for a in range(n):
        for b in range(n):
            F[a, b] = fp(lam[a]) if a == b else (f(lam[a]) - f(lam[b])) / (lam[a] - lam[b])
    Uh = U.transpose_conj()
    D = []
    for g in G:
        Mk = Uh * (g / 2) * U
        for a in range(n):
            for b in range(n): Mk[a, b] *= F[a, b]
        D.append(U * Mk * Uh)
    m = len(G)
    gm = mp.matrix(m, m)
    for i in range(m):
        for j in range(i, m):
            v = mp.re(sum((D[i] * D[j])[a, a] for a in range(n)))
            gm[i, j] = v; gm[j, i] = v
    return gm, lam

def scalar_curvature(n, x0, f, fp, h=mp.mpf('1e-15')):
    G = herm_basis(n); m = len(G)
    def g_at(x): return metric(n, G, x, f, fp)[0]
    def gamma_at(x):
        g = g_at(x); gi = g ** -1
        dg = []
        for l in range(m):
            xp = list(x); xm = list(x); xp[l] += h; xm[l] -= h
            dg.append((g_at(xp) - g_at(xm)) / (2 * h))
        Gam = [[[sum(gi[k, l] * (dg[j][l, i] + dg[i][l, j] - dg[l][i, j]) for l in range(m)) / 2
                 for j in range(m)] for i in range(m)] for k in range(m)]
        return Gam, gi
    Gam, gi = gamma_at(x0)
    dGam = []
    for l in range(m):
        xp = list(x0); xm = list(x0); xp[l] += h; xm[l] -= h
        Gp = gamma_at(xp)[0]; Gm = gamma_at(xm)[0]
        dGam.append([[[(Gp[k][i][j] - Gm[k][i][j]) / (2 * h) for j in range(m)] for i in range(m)] for k in range(m)])
    R = 0
    for i in range(m):
        for j in range(m):
            Ric = sum(dGam[k][k][i][j] - dGam[j][k][i][k]
                      + sum(Gam[k][k][l] * Gam[l][i][j] - Gam[k][j][l] * Gam[l][i][k] for l in range(m))
                      for k in range(m))
            R += gi[i, j] * Ric
    lam = metric(n, G, x0, f, fp)[1]
    return R, lam

f3 = lambda t: 3 * mp.cbrt(t); f3p = lambda t: mp.cbrt(t) / t
f2 = lambda t: 2 * mp.sqrt(t); f2p = lambda t: 1 / mp.sqrt(t)
for n, x0 in [(2, [mp.mpf('0.2'), mp.mpf('0.1'), mp.mpf('0.3')]),
              (3, [mp.mpf(v) for v in ('0.05', '-0.03', '0.08', '0.02', '-0.06', '0.04', '0.11', '-0.07')])]:
    R, lam = scalar_curvature(n, x0, f3, f3p)
    P = pred_p3(lam)
    print(f'(1) n={n} intrinsic Scal_3 = {mp.nstr(R, 25)}  N/(9 s4^2) = {mp.nstr(P, 25)}')
    report(f'(1) n={n} intrinsic curvature equals the formula:', abs(R - P) < mp.mpf('1e-20') * max(1, abs(P)))
    R2, _ = scalar_curvature(n, x0, f2, f2p)
    e2 = mp.mpf((n*n - 1) * (n*n - 2)) / 4
    print(f'(1) n={n} control p=2: {mp.nstr(R2, 25)}  expected {e2}')
    report(f'(1) n={n} control p=2 equals (n^2-1)(n^2-2)/4:', abs(R2 - e2) < mp.mpf('1e-20') * e2)

# ---------------- (2) Gauss equation with an explicit shape-operator matrix ----------------
mp.mp.dps = 40
random.seed(7)
def orth_herm_basis(n):
    B = []
    for i in range(n):
        M = mp.zeros(n); M[i, i] = 1; B.append(M)
    for i in range(n):
        for j in range(i + 1, n):
            M = mp.zeros(n); M[i, j] = 1 / mp.sqrt(2); M[j, i] = 1 / mp.sqrt(2); B.append(M)
            M = mp.zeros(n); M[i, j] = 1j / mp.sqrt(2); M[j, i] = -1j / mp.sqrt(2); B.append(M)
    return B
def ip(A, B, n): return mp.re(sum((A * B)[a, a] for a in range(n)))
def rand_unitary(n):
    Z = mp.matrix(n, n)
    for a in range(n):
        for b in range(n): Z[a, b] = mp.mpf(random.gauss(0, 1)) + 1j * mp.mpf(random.gauss(0, 1))
    Q, _ = mp.qr(Z)
    return Q
for n in range(2, 6):
    for trial in range(3):
        ys = [mp.mpf(random.uniform(0.2, 1.5)) for _ in range(n)]
        sc = mp.cbrt(sum(t**3 for t in ys)); ys = [t / sc for t in ys]
        U = rand_unitary(n); D = mp.diag(ys); Y = U * D * U.transpose_conj()
        B = orth_herm_basis(n); m = n * n
        K = mp.matrix(m, m)
        for a in range(m):
            for b in range(m): K[a, b] = ip(B[a], Y * B[b] + B[b] * Y, n)
        s4 = sum(t**4 for t in ys)
        Y2 = Y * Y
        nu = mp.matrix([ip(B[a], Y2, n) / mp.sqrt(s4) for a in range(m)])
        P = mp.eye(m) - nu * nu.T
        S = P * K * P / mp.sqrt(s4)
        trS = sum(S[a, a] for a in range(m)); trS2 = sum((S * S)[a, a] for a in range(m))
        gauss = trS**2 - trS2
        pred = N_of(ys, n) / s4**2
        if trial == 0: print(f'(2) n={n} (trS)^2 - tr S^2 = {mp.nstr(gauss, 20)}  N/s4^2 = {mp.nstr(pred, 20)}')
        report(f'(2) n={n} trial {trial} Gauss evaluation equals N/s4^2:', abs(gauss - pred) < mp.mpf('1e-25') * max(1, abs(pred)))

# ---------------- (3) symbolic algebra ----------------
x, y, a, b, c, d, e, n = sp.symbols('x y a b c d e n')
s1, s2, s4, s5, s6 = x + y + a, x**2 + y**2 + b, x**4 + y**4 + c, x**5 + y**5 + d, x**6 + y**6 + e
Nsym = ((4*n**2 - 2)*s1**2 - 2*n*s2)*s4 - 8*n*s1*s5 + 8*s6
num = x**2*(sp.diff(Nsym, y)*s4 - 8*Nsym*y**3) - y**2*(sp.diff(Nsym, x)*s4 - 8*Nsym*x**3)
Qsym, rem = sp.div(sp.expand(num), x - y, x)
report('(3) numerator divisible by x - y:', sp.expand(rem) == 0)
Qsym = sp.expand(Qsym / 4)
for nn in (2, 3, 4):
    ys = sp.symbols(f'y1:{nn + 1}', positive=True)
    S = lambda k: sum(t**k for t in ys)
    Nn = ((4*nn**2 - 2)*S(1)**2 - 2*nn*S(2))*S(4) - 8*nn*S(1)*S(5) + 8*S(6)
    R = Nn / (9 * S(4)**2)
    i, j = 0, 1
    lhs = sp.diff(R, ys[j]) / (3 * ys[j]**2) - sp.diff(R, ys[i]) / (3 * ys[i]**2)
    zs = ys[2:]
    sub = {x: ys[i], y: ys[j], a: sum(zs), b: sum(t**2 for t in zs), c: sum(t**4 for t in zs),
           d: sum(t**5 for t in zs), e: sum(t**6 for t in zs), n: nn}
    rhs = 4 * (ys[i] - ys[j]) * Qsym.subs(sub) / (27 * ys[i]**2 * ys[j]**2 * S(4)**3)
    report(f'(3) chain-rule identity n={nn}:', sp.simplify(sp.together(lhs - rhs)) == 0)
def H(m, *cs):
    out = sum(cs[k] * (x**(m - k) * y**k + x**k * y**(m - k)) for k in range(m) if 2 * k < m)
    if m % 2 == 0: out += cs[m // 2] * x**(m // 2) * y**(m // 2)
    return out
u, v, w = a**2 - b, a*c - d, e - 2*a*d + a**2*c
Lsym = sp.expand(Qsym.subs({b: a**2, d: a*c, e: a**2*c}, simultaneous=True))
L_disp = ((4*n**2 - 18*n + 14)*a**2*c*x**2*y**2 + (4*n**2 - 2*n - 2)*a**2*x**2*y**2*(x**4 + y**4)
          + (2*n**2 - 2*n - 1)*a*c**2*(x + y) + a*c*H(5, 4*n**2 - 2*n - 2, 4*n**2 - 2*n - 2, 8*n**2 - 6*n - 4)
          + a*H(9, 2*n**2 - 1, 2*n**2 - 1, 8*n**2 - 6*n - 4, 8*n**2 + 10*n - 4, 4*n**2 - 2)
          + c**2*H(2, 2*n**2 - 1, 4*n**2 - n - 2)
          + c*H(6, 4*n**2 - 2*n - 2, 8*n**2 - 4*n - 4, 8*n**2 + 8*n - 16, 8*n**2 + 20*n - 16)
          + H(10, 2*n**2 - 2*n - 1, 4*n**2 - 3*n - 2, 6*n**2 - 8*n + 1, 8*n**2 + 4*n - 16, 8*n**2 + 6*n - 16, 8*n**2 - 6*n - 4))
report('(3) L equals the displayed closed form:', sp.expand(Lsym - L_disp) == 0)
deficit = (2*n*x**2*y**2*(c + x**4 + y**4)*u + 2*(8*a*(n - 2)*x**2*y**2 + n*c*(x + y) + n*H(5, 1, 1, 8))*v
           + 16*x**2*y**2*w)
report('(3) moment-deficit identity Q = L + deficits (symbolic n):', sp.expand(Qsym - Lsym - deficit) == 0)
k, z = sp.symbols('k z')
Pk = sp.Poly(sp.expand(L_disp.subs(n, k + 4)), k, a, c, x, y)
cs = Pk.coeffs()
print('(3) n=k+4: monomials', len(cs), 'min coefficient', min(cs))
report('(3) n>=4: every coefficient of L(n=k+4) is positive:', all(t > 0 for t in cs))
P2 = sp.Poly(sp.expand(L_disp.subs({n: 2, a: 0, c: 0})), x, y)
report('(3) n=2: L = H_10(3,8,9,24,28,16):', sp.expand(P2.as_expr() - H(10, 3, 8, 9, 24, 28, 16)) == 0)
print('(3) n=2: monomials', len(P2.coeffs()))
report('(3) n=2: coefficients positive:', all(t > 0 for t in P2.coeffs()))
L3 = sp.expand(L_disp.subs({n: 3, a: z, c: z**4}))
neg = [t for t in sp.Poly(L3, x, y, z).terms() if t[1] < 0]
print('(3) n=3: negative monomials of L', neg)
P3 = sp.Poly(sp.expand(L3 - y**2*z**4*(2*x**2 - z**2)**2), x, y, z)
print('(3) n=3: monomials after subtracting the square', len(P3.coeffs()))
report('(3) n=3: remainder coefficients positive:', all(t > 0 for t in P3.coeffs()))
report('(3) n=3: x^10 present in the remainder:', P3.coeff_monomial(x**10) > 0)

# ---------------- (4) strict Schur monotonicity under T-transforms ----------------
random.seed(11)
def Rf(lams):
    ys = [t ** (1.0 / 3) for t in lams]; nn = len(ys)
    return N_of(ys, nn) / (9 * sum(t**4 for t in ys)**2)
viol = 0; trials = 0
for _ in range(20000):
    nn = random.randint(2, 7)
    lam = [random.random()**3 + 1e-3 for _ in range(nn)]; s = sum(lam); lam = [t / s for t in lam]
    i, j = random.sample(range(nn), 2); t = random.uniform(0.01, 0.99)
    mu = list(lam); mu[i] = t*lam[i] + (1 - t)*lam[j]; mu[j] = (1 - t)*lam[i] + t*lam[j]
    if abs(lam[i] - lam[j]) > 1e-6:
        trials += 1
        if not Rf(mu) > Rf(lam): viol += 1
print('(4) nontrivial T-transforms', trials, 'violations', viol)
report('(4) no violation:', viol == 0)
print('bad=', bad)
sys.exit(1 if bad else 0)
