"""Checks for Bera-Sarbicki-Chruscinski, arXiv:2309.09621v1, Conjecture 2 Eq. (13): for GCD(n,k)=3 and |alpha|=|beta|, the largest
lambda with tau_{n,k} - lambda H_v positive would be n-k+(n-3)/(n-2k/3). Exact rational counterexample at (n,k)=(6,3), alpha=beta=1/sqrt2;
the residue-class compression x = 1_{3Z} for every (n,k) with GCD 3, n <= 30; numerical lambda_max at equator phases and at the pole for (6,3)."""
from math import gcd
import sympy as sp, numpy as np

def tau(X, n, k):
    E = lambda Y: sp.diag(*[Y[i, i] for i in range(n)])
    S = sp.zeros(n)
    for i in range(n):
        S[(i + 1) % n, i] = 1
    out = (n - k) * E(X) - X
    for a in range(1, k + 1):
        Sa = S ** a
        out += E(Sa * X * Sa.H)
    return out

ok = True
n, k = 6, 3
z = sp.exp(2 * sp.pi * sp.I / 3)
v = (sp.Matrix([z ** i for i in range(n)]) + sp.Matrix([z ** (2 * i) for i in range(n)])) / sp.sqrt(n) / sp.sqrt(2)
v = sp.Matrix([sp.nsimplify(sp.simplify(sp.expand_complex(x))) for x in v])
P = v * v.H
lam = sp.Rational(n - k) + sp.Rational(n - 3) / (n - sp.Rational(2, 3) * k)
q = sp.Matrix([1, 0, 0, 1, 0, 0]); R = (sp.eye(6) + 27 * q * q.T) / 60
PR = sp.simplify(tau(R, n, k) - lam * P.multiply_elementwise(R))
val = sp.simplify((q.T * PR * q)[0])
c1 = lam == sp.Rational(15, 4) and R.trace() == 1 and min(R.eigenvals()) > 0 and val == sp.Rational(-13, 40)
print("exact counterexample (6,3): lambda* = %s, q* Phi(R) q = %s, R positive definite with trace 1:" % (lam, val), c1)
ok &= c1

fam = True; cases = 0
for N in range(6, 31, 3):
    for K in range(1, N):
        if gcd(N, K) != 3:
            continue
        cases += 1
        vv = (np.array([np.exp(2j*np.pi*i/3) for i in range(N)]) + np.array([np.exp(4j*np.pi*i/3) for i in range(N)])) / np.sqrt(2*N)
        Pn = np.outer(vv, vv.conj()); x = np.array([1.0 if i % 3 == 0 else 0.0 for i in range(N)]); X = np.outer(x, x)
        Sn = np.roll(np.eye(N), 1, axis=0); out = (N - K) * np.diag(np.diag(X)) - X
        for a in range(1, K + 1):
            Sa = np.linalg.matrix_power(Sn, a); out = out + np.diag(np.diag(Sa @ X @ Sa.T))
        lstar = N - K + (N - 3) / (N - 2 * K / 3)
        Y = out - lstar * Pn * X
        fam &= np.linalg.eigvalsh(Y).min() < -1e-9 and lstar > N - K
print("family: x = 1_{3Z} certifies non-positivity at the conjectured value for all %d pairs (n,k), GCD 3, n <= 30:" % cases, fam)
ok &= fam

rng = np.random.default_rng(7)
def lammax(Pn, N=6, K=3):
    Sn = np.roll(np.eye(N), 1, axis=0); Sp = [np.linalg.matrix_power(Sn, a) for a in range(1, K + 1)]
    def taun(X):
        out = (N - K) * np.diag(np.diag(X)) - X
        for Sa in Sp:
            out = out + np.diag(np.diag(Sa @ X @ Sa.conj().T))
        return out
    def f(x, l):
        x = x / np.linalg.norm(x); X = np.outer(x, x.conj()); return np.linalg.eigvalsh(taun(X) - l * Pn * X).min()
    def worst(l):
        best = 1e9
        for _ in range(30):
            x = rng.normal(size=N) + 1j * rng.normal(size=N); fx = f(x, l); step = 0.5
            for _ in range(1200):
                y = x + step * (rng.normal(size=N) + 1j * rng.normal(size=N)); fy = f(y, l)
                if fy < fx: x, fx = y, fy
                else: step *= 0.995
            best = min(best, fx)
        return best
    lo, hi = 2.5, 4.5
    for _ in range(11):
        mid = (lo + hi) / 2
        if worst(mid) >= -1e-7: lo = mid
        else: hi = mid
    return lo
v1 = np.array([np.exp(2j*np.pi*i/3) for i in range(6)]) / np.sqrt(6); v2 = np.array([np.exp(4j*np.pi*i/3) for i in range(6)]) / np.sqrt(6)
for phi in (0.0, np.pi / 3, np.pi):
    w = (v1 + np.exp(1j*phi) * v2) / np.sqrt(2)
    print("numerical lambda_max (6,3), |alpha|=|beta|, phase %.4f: %.4f" % (phi, lammax(np.outer(w, w.conj()))))
print("numerical lambda_max (6,3), alpha=1, beta=0: %.4f" % lammax(np.outer(v1, v1.conj())))
print("RESULT", "PASS" if ok else "FAIL")
