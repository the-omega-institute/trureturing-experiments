#!/usr/bin/env python3
"""Numerical check of the proposed answer to Zhan-Gupta-Hillery (ACM TQC 5(2), 2024; arXiv:2306.17401v6) Conjecture 1.

Setting: n qubit sensors, U has eigenvalues e^{-i theta}, e^{+i theta}; final states phi_i = U_i psi, equal priors 1/n.
Claim: for theta in (0,T] (T = arccos(-(q-1)/q)/2, q = ceil(n/2)), the Dicke state D_l (l = floor(n/2)) is optimal,
with optimal success d(c*)^2, d(c) = (sqrt(1+(n-1)c) + (n-1) sqrt(1-c))/n, c* = 1 - 4 floor(n^2/4) sin^2(theta)/(n(n-1)).

Checks (all exact formulas against an independent iterative optimum):
 1. Gram matrix of U_i psi depends only on |psi_x|^2 (random phases).
 2. Optimal success of the Dicke ensemble (iterative maximum-success POVM) equals d(c*)^2.
 3. For random inputs psi, the iterative optimum never exceeds d(c*)^2.
 4. The dual certificate Y = (d/n)H satisfies Y >= H e_i e_i^* H / n for constant-overlap Gram H^2.
"""
import itertools, math, sys
import numpy as np

rng = np.random.default_rng(20261014)

def finals(n, theta, psi):
    xs = np.array(list(itertools.product([0, 1], repeat=n)))
    out = []
    for i in range(n):
        ph = np.exp(1j * theta * (2 * xs[:, i] - 1))  # bit 0 -> e^{-i theta}, bit 1 -> e^{+i theta}
        out.append(ph * psi)
    return np.array(out)  # n x 2^n

def gram(V):
    return V.conj() @ V.T  # G_ij = <v_i, v_j>

def opt_success(G, iters=4000):
    """Maximum average success for equal-prior pure states with Gram G (Jezek-Rehacek-Fiurasek iteration in span)."""
    n = G.shape[0]
    w, Q = np.linalg.eigh(G)
    w = np.clip(w, 0, None)
    S = Q @ np.diag(np.sqrt(w)) @ Q.conj().T  # columns are vectors with Gram G (G = S^* S since S Hermitian)
    rho = [np.outer(S[:, i], S[:, i].conj()) / n for i in range(n)]
    Pi = [np.eye(n) / n for _ in range(n)]
    for _ in range(iters):
        R = sum(rho[i] @ Pi[i] @ rho[i] for i in range(n))
        wr, Qr = np.linalg.eigh(R)
        wr = np.where(wr > 1e-14, wr, np.inf)
        Rm = Qr @ np.diag(1 / np.sqrt(wr)) @ Qr.conj().T
        Pi = [Rm @ rho[i] @ Pi[i] @ rho[i] @ Rm for i in range(n)]
        # complete on the kernel of R so that sum Pi = I (does not change success)
    succ = sum(np.real(np.trace(rho[i] @ Pi[i])) for i in range(n))
    return succ

def dval(n, c):
    return ((math.sqrt(1 + (n - 1) * c) + (n - 1) * math.sqrt(1 - c)) / n) ** 2

def dicke(n, l):
    xs = list(itertools.product([0, 1], repeat=n))
    v = np.array([1.0 if sum(x) == l else 0.0 for x in xs], dtype=complex)
    return v / np.linalg.norm(v)

ok = True
for n in (3, 4, 5):
    q = math.ceil(n / 2)
    T = math.acos(-(q - 1) / q) / 2
    for frac in (0.25, 0.5, 0.8, 1.0):
        theta = frac * T
        cstar = 1 - 4 * (n * n // 4) * math.sin(theta) ** 2 / (n * (n - 1))
        best = dval(n, cstar)
        # 1. phase independence of the Gram matrix
        p = rng.dirichlet(np.ones(2 ** n))
        ph1 = np.exp(2j * np.pi * rng.random(2 ** n)); ph2 = np.exp(2j * np.pi * rng.random(2 ** n))
        G1 = gram(finals(n, theta, np.sqrt(p) * ph1)); G2 = gram(finals(n, theta, np.sqrt(p) * ph2))
        e1 = np.max(np.abs(G1 - G2))
        # 2. Dicke optimum
        GD = gram(finals(n, theta, dicke(n, n // 2)))
        sD = opt_success(GD)
        # 3. random inputs
        worst = 0.0
        for _ in range(12):
            alpha = rng.choice([0.2, 1.0, 5.0])
            pr = rng.dirichlet(alpha * np.ones(2 ** n))
            psi = np.sqrt(pr) * np.exp(2j * np.pi * rng.random(2 ** n))
            s = opt_success(gram(finals(n, theta, psi)), iters=1500)
            worst = max(worst, s)
        # 4. dual certificate for constant overlap c* (c* >= 0 on the range)
        J = np.ones((n, n)); P = J / n; Qm = np.eye(n) - P
        H = math.sqrt(1 + (n - 1) * cstar) * P + math.sqrt(1 - cstar) * Qm
        d = math.sqrt(best)
        mineig = min(np.linalg.eigvalsh((d / n) * H - np.outer(H[:, i], H[:, i]) / n).min() for i in range(n))
        line_ok = e1 < 1e-12 and abs(sD - best) < 1e-6 and worst <= best + 1e-9 and mineig > -1e-12 and cstar >= -1e-12
        ok &= line_ok
        print(f"n={n} theta/T={frac:.2f} c*={cstar:.6f} d(c*)^2={best:.9f} dicke_opt={sD:.9f} "
              f"max_random_opt={worst:.9f} gram_phase_err={e1:.1e} cert_min_eig={mineig:.1e} {'OK' if line_ok else 'FAIL'}")
print("ALL_OK" if ok else "SOME_FAIL")
sys.exit(0 if ok else 1)
