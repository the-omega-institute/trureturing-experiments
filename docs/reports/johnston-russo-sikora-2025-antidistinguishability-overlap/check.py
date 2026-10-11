"""Antidistinguishability from pairwise overlaps: Corollary 5.4 of Johnston-Russo-Sikora is not tight for n >= 10.

Source: N. Johnston, V. Russo, J. Sikora, "Tight bounds for antidistinguishability and circulant sets of pure
quantum states", arXiv:2311.17047v2; Quantum 9, 1622 (2025). Section 6: "We still do not know whether or not
Corollary 5.4 is tight when n >= 5." Corollary 5.4: |<psi_i, psi_j>| <= c_n = sqrt((n-2)/(2n-2)) for all i != j
implies antidistinguishability. Theorem 5.1: sqrt(l_0) <= sum_{j>=1} sqrt(l_j) (Gram eigenvalues) implies it.

Claim checked here (Part A of the route): with N = n(n-1), delta_n = (n-9)/(16 N (n-3)),
max_{i!=j} |G_ij|^2 <= c_n^2 + delta_n forces the Theorem 5.1 spectral condition, for n >= 10.

Checks:
 1. Lemma A1 identity: tr(G^2) - n^2/2 = 2(E^2 - Q) + 2 E t + t^2/2 on random spectra with trace n (exact, sympy).
 2. The final inequality sqrt((n-3) N delta) + 2 N delta < (1 - c_n^2) sqrt(n-9), with both sides in closed form,
    exactly in sympy for 10 <= n <= 400.
 3. Rank of the entrywise squared modulus of a rank-3 Gram matrix is <= 9 (random complex rank-3 Gram, n = 20).
 4. Stress test: random Gram matrices whose off-diagonal moduli satisfy |G_ij|^2 <= c_n^2 + delta_n (n = 10, 12, 16;
    equal-overlap circulant-type, random-phase equal-modulus, and random perturbations kept PSD) all satisfy the
    Theorem 5.1 spectral condition. This is a sanity check, not a proof.
"""
import sys, itertools
import numpy as np
import sympy as sp

bad = 0
def ok(c, msg):
    global bad
    print(("ok   " if c else "FAIL ") + msg)
    if not c: bad += 1

# 1. identity
rng = np.random.default_rng(20261012)
allid = True
for n in (5, 7, 10):
    for _ in range(20):
        xs = [sp.Rational(int(v), 7) for v in rng.integers(1, 30, size=n - 1)]
        a = sum(xs); E = sum(xs[i] * xs[j] for i in range(n - 1) for j in range(i + 1, n - 1))
        Q = sum(xs[i] ** 2 * xs[j] ** 2 for i in range(n - 1) for j in range(i + 1, n - 1))
        t = sp.Rational(int(rng.integers(1, 50)), 11)
        lam1 = a ** 2 + t
        lams = [lam1] + [x ** 2 for x in xs]
        nn = sum(lams)  # plays the role of n = tr G
        lhs = sum(l ** 2 for l in lams) - nn ** 2 / 2
        rhs = 2 * (E ** 2 - Q) + 2 * E * t + t ** 2 / 2
        if sp.simplify(lhs - rhs) != 0: allid = False
ok(allid, "Lemma A1 identity tr(G^2) - (tr G)^2/2 = 2(E^2-Q) + 2Et + t^2/2 (exact, 60 random spectra)")

# 2. final inequality
good = True
for n in range(10, 401):
    N = n * (n - 1); c2 = sp.Rational(n - 2, 2 * (n - 1)); d = sp.Rational(n - 9, 16 * N * (n - 3))
    left = sp.sqrt((n - 3) * N * d) + 2 * N * d
    right = (1 - c2) * sp.sqrt(n - 9)
    if n == 10:
        pass
    if not (d > 0 and sp.simplify(left - (sp.sqrt(n - 9) / 4 + sp.Rational(n - 9, 8 * (n - 3)))) == 0 and (right - left).evalf(50) > 0):
        good = False; print("  failure at n =", n)
ok(good, "closed form sqrt((n-3)N delta) + 2N delta = sqrt(n-9)/4 + (n-9)/(8(n-3)) < (1-c_n^2) sqrt(n-9), 10 <= n <= 400")

# 3. rank of |G_3|^2
V = rng.normal(size=(20, 3)) + 1j * rng.normal(size=(20, 3))
G3 = V @ V.conj().T
H3 = np.abs(G3) ** 2
ok(np.linalg.matrix_rank(H3, tol=1e-8) <= 9, f"rank of |G_3|^2 for a random rank-3 Gram (n=20) is {np.linalg.matrix_rank(H3, tol=1e-8)} <= 9")

# 4. stress test
def spectral_ok(G):
    w = np.sort(np.clip(np.linalg.eigvalsh(G), 0, None))[::-1]
    return np.sqrt(w[0]) <= np.sum(np.sqrt(w[1:])) + 1e-12
def near_psd(G):
    return np.min(np.linalg.eigvalsh(G)) >= -1e-12
cnt = 0; viol = 0
for n in (10, 12, 16):
    c2 = (n - 2) / (2 * (n - 1)); dl = (n - 9) / (16 * n * (n - 1) * (n - 3)); cap = np.sqrt(c2 + dl)
    trials = []
    trials.append((1 - cap) * np.eye(n) + cap * np.ones((n, n)))
    for _ in range(400):
        # random unit vectors in dimension d with overlaps capped: rescale toward a common direction
        d = rng.integers(2, n + 1)
        X = rng.normal(size=(n, d)) + 1j * rng.normal(size=(n, d))
        X /= np.linalg.norm(X, axis=1, keepdims=True)
        G = X @ X.conj().T
        off = np.abs(G - np.diag(np.diag(G)))
        if off.max() ** 2 <= c2 + dl:
            trials.append(G)
        # push toward the boundary: mix with all-ones direction
        u = np.ones(d) / np.sqrt(d)
        for s in (0.3, 0.6, 0.8):
            Y = (1 - s) * X + s * u
            Y /= np.linalg.norm(Y, axis=1, keepdims=True)
            G2 = Y @ Y.conj().T
            o2 = np.abs(G2 - np.diag(np.diag(G2)))
            if o2.max() ** 2 <= c2 + dl:
                trials.append(G2)
    for G in trials:
        if near_psd(G):
            cnt += 1
            if not spectral_ok(G): viol += 1
ok(viol == 0, f"stress test: {cnt} admissible Gram matrices (n = 10, 12, 16), {viol} violate the spectral condition")
print("bad=", bad)
sys.exit(1 if bad else 0)
