"""Cha–Lee arXiv:2602.17748v3, Problem 1 (last question): strict improvement alpha* < 1/sqrt(2).
Proposed proof: ||Gamma X||_1 <= (70/99) d ||X||_1 for Hermitian X on C^d (x) C^d with Tr_A X = 0, Gamma the partial
transpose on A; hence ||Theta o Psi||_dia <= (70/99) d ||Psi||_dia for trace-annihilating Hermiticity-preserving Psi.
1. Lemma 1: the polynomial identity q(r) - r = (2-r)(r-1)^2(r+4)/18 and the rewriting of q (sympy); random tests of
   sum a_i a_j |z_i - z_j| <= (4/3)(sum a_i)^2.
2. Lemma 2: ||Gamma(P_u - P_v)||_1 <= 4d/3 for random equal-marginal pairs v = (U (x) I) u.
3. Matrix theorem: random and adversarially optimized ratios ||Gamma X||_1 / (d ||X||_1) for d = 2, 3 stay below 70/99.
4. Exact rational comparisons of Lemma 4 (fractions) and the source's cyclic lower bound cot(pi/2d)/d < 70/99.
"""
import math, sys
from fractions import Fraction as F
import numpy as np
import sympy as sp
from scipy.optimize import minimize

rng = np.random.default_rng(20261017)


def ptA(X, d):  # partial transpose on the first factor of C^d (x) C^d
    return X.reshape(d, d, d, d).transpose(2, 1, 0, 3).reshape(d * d, d * d)


def trA(X, d):
    return np.einsum('ijik->jk', X.reshape(d, d, d, d))


def tn(X):
    return np.abs(np.linalg.eigvalsh((X + X.conj().T) / 2)).sum()


def tn_general(X):
    return np.linalg.svd(X, compute_uv=False).sum()


def main():
    fails = []
    r, x, c2 = sp.symbols('r x c2', real=True)
    q = sp.Rational(4, 9) + sp.Rational(11, 18) * r ** 2 - r ** 4 / 18
    if sp.expand(q - r - (2 - r) * (r - 1) ** 2 * (r + 4) / 18) != 0: fails.append("lemma1-identity")
    # with r^2 = 2 - 2x and Re(w^2) = 2x^2 - 1
    q_x = sp.expand(q.subs(r ** 2, 2 - 2 * x).subs(r ** 4, (2 - 2 * x) ** 2))
    q_x = sp.expand(sp.Rational(4, 9) + sp.Rational(11, 18) * (2 - 2 * x) - (2 - 2 * x) ** 2 / 18)
    if sp.expand(q_x - (sp.Rational(4, 3) - sp.Rational(7, 9) * x - sp.Rational(1, 9) * (2 * x ** 2 - 1))) != 0: fails.append("lemma1-rewrite")
    worst1 = 0
    for _ in range(20000):
        n = rng.integers(1, 9); a = rng.random(n); z = np.exp(1j * rng.uniform(0, 2 * np.pi, n))
        if rng.random() < 0.3: z = np.exp(2j * np.pi * np.arange(n) / n)
        lhs = sum(a[i] * a[j] * abs(z[i] - z[j]) for i in range(n) for j in range(n))
        worst1 = max(worst1, lhs / (a.sum() ** 2))
    print(f"Lemma 1: identity ok = {'lemma1-identity' not in fails and 'lemma1-rewrite' not in fails}; max ratio sum a a |z-z| / (sum a)^2 = {worst1:.5f} (bound 4/3 = 1.33333; 4/pi = {4/math.pi:.5f})")
    if worst1 > 4 / 3 + 1e-12: fails.append("lemma1")
    worst2 = 0
    for _ in range(3000):
        d = int(rng.integers(2, 6))
        u = rng.normal(size=d * d) + 1j * rng.normal(size=d * d); u /= np.linalg.norm(u)
        Z = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)); U, _ = np.linalg.qr(Z)
        v = np.kron(U, np.eye(d)) @ u
        if np.abs(trA(np.outer(u, u.conj()), d) - trA(np.outer(v, v.conj()), d)).max() > 1e-12: fails.append("lemma2-marginal"); break
        val = tn_general(ptA(np.outer(u, u.conj()) - np.outer(v, v.conj()), d))
        worst2 = max(worst2, val / d)
    print(f"Lemma 2: max ||Gamma(P_u - P_v)||_1 / d over 3000 random equal-marginal pairs = {worst2:.5f} (bound 4/3)")
    if worst2 > 4 / 3 + 1e-9: fails.append("lemma2")

    def ratio_from(params, d):
        n = d * d
        Y = (params[:n * n] + 1j * params[n * n:]).reshape(n, n); Y = Y + Y.conj().T
        X = Y - np.kron(np.eye(d), trA(Y, d)) / d      # Tr_A X = 0
        nx = tn(X)
        return 0.0 if nx < 1e-12 else tn_general(ptA(X, d)) / (d * nx)

    for d in (2, 3):
        best = 0
        for _ in range(200):
            p = rng.normal(size=2 * d ** 4); best = max(best, ratio_from(p, d))
        for _ in range(12 if d == 2 else 6):
            res = minimize(lambda p: -ratio_from(p, d), rng.normal(size=2 * d ** 4), method='Nelder-Mead',
                           options={'maxiter': 20000 if d == 2 else 30000, 'xatol': 1e-9, 'fatol': 1e-12})
            best = max(best, -res.fun)
        print(f"matrix theorem, d={d}: max ratio ||Gamma X||_1 / (d ||X||_1) found = {best:.6f} (70/99 = {70/99:.6f}, 1/sqrt2 = {1/math.sqrt(2):.6f})")
        if best > 70 / 99: fails.append(("matrix", d, best))
    # control: Psi = id - Ad_U for the cyclic shift U; X = (Psi (x) id)(maximally entangled state) and the best input
    for d in (2, 3, 4, 6, 8):
        U = np.roll(np.eye(d), 1, axis=0)
        best_in = 0
        for trial in range(60 if d <= 4 else 20):
            if trial == 0:
                psi = np.eye(d).reshape(d * d) / math.sqrt(d)
            else:
                psi = rng.normal(size=d * d) + 1j * rng.normal(size=d * d); psi /= np.linalg.norm(psi)
            P = np.outer(psi, psi.conj()); UU = np.kron(U, np.eye(d))
            X = P - UU @ P @ UU.conj().T
            if tn(X) < 1e-12: continue
            best_in = max(best_in, tn_general(ptA(X, d)) / (d * tn(X)))
        cyc = 1 / math.tan(math.pi / (2 * d)) / d
        print(f"cyclic control d={d}: ratio at max-entangled/random inputs = {best_in:.6f}; source cot(pi/2d)/d = {cyc:.6f}")
        if best_in > 70 / 99: fails.append(("cyclic", d))
    ok4a = F(8191) * 9801 < F(4900) * 16384
    ok4b = F(2, 3) + 2 * F(1, 4096) + 2 * F(1, 64) == F(4291, 6144) and F(4291, 6144) < F(70, 99)
    ok70 = 2 * F(70) ** 2 < F(99) ** 2
    print(f"exact comparisons: case 4a {ok4a}, case 4b {ok4b}, (70/99)^2 < 1/2 {ok70}")
    if not (ok4a and ok4b and ok70): fails.append("rational")
    lb = max(1 / math.tan(math.pi / (2 * d)) / d for d in range(2, 5000))
    print(f"source lower bound sup_d cot(pi/2d)/d = {lb:.6f} (-> 2/pi = {2/math.pi:.6f}) < 70/99: {lb < 70/99}")
    if lb >= 70 / 99: fails.append("lower-bound")
    print("RESULT", "PASS" if not fails else f"FAIL {fails}")
    sys.exit(0 if not fails else 1)


if __name__ == "__main__":
    main()
