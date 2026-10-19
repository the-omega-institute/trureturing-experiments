#!/usr/bin/env python3
"""Louvet-Serrano-Bastin-Martin (arXiv:2411.16461v2) Conjecture 1: numerical checks.

For N sites of dimension d, P is the projector onto Sym^N(C^d), D = C(N+d-1, d-1).
Checks, for every subset S of sites with |S| = k:
  (1) lambda_min of (P/D)^{T_S} restricted to Sym^k x Sym^{N-k} equals 1/(D C(N,k)),
      and P^{T_S} - Pi_S x Pi_Sc / C(N,k) is PSD;
  (2) lambda_min((|psi><psi|)^{T_S}) >= -1/2 for random unit psi;
  (3) rho = (p/D) P + (1-p)|psi><psi| has PSD partial transpose for every S at p = p_min,
      for random symmetric psi and for the GHZ-type psi; and it is NPT below p_min for the GHZ-type psi.
Exit 0 and print RESULT PASS when all checks hold.
"""
import itertools, math, sys
import numpy as np

rng = np.random.default_rng(20261019)
TOL = 1e-9

def perm_op(N, d, sigma):
    n = d ** N
    M = np.zeros((n, n))
    for idx, x in enumerate(itertools.product(range(d), repeat=N)):
        y = tuple(x[sigma[i]] for i in range(N))
        jdx = int(np.ravel_multi_index(y, (d,) * N)) if N else 0
        M[jdx, idx] = 1
    return M

def sym_proj(N, d):
    n = d ** N
    P = np.zeros((n, n))
    for s in itertools.permutations(range(N)):
        P += perm_op(N, d, s)
    return P / math.factorial(N)

def pt(rho, N, d, S):
    t = rho.reshape((d,) * (2 * N))
    axes = list(range(2 * N))
    for i in S:
        axes[i], axes[N + i] = axes[N + i], axes[i]
    return t.transpose(axes).reshape(d ** N, d ** N)

def sub_sym(N, d, S):
    """Projector onto Sym(sites in S) x Sym(sites outside S)."""
    Sc = [i for i in range(N) if i not in S]
    P = np.zeros((d ** N, d ** N))
    cnt = 0
    for a in itertools.permutations(S):
        for b in itertools.permutations(Sc):
            sigma = list(range(N))
            for i, j in zip(S, a): sigma[i] = j
            for i, j in zip(Sc, b): sigma[i] = j
            P += perm_op(N, d, sigma); cnt += 1
    return P / cnt

ok = True
for N, d in [(2, 2), (2, 3), (3, 2), (3, 3), (4, 2), (2, 4), (4, 3), (5, 2), (3, 4)]:
    P = sym_proj(N, d)
    D = math.comb(N + d - 1, d - 1)
    assert abs(np.trace(P) - D) < TOL
    C = math.comb(N, N // 2)
    pmin = 1 / (1 + 2 / (D * C))
    w, V = np.linalg.eigh(P)
    symbasis = V[:, w > 0.5]
    for k in range(0, N + 1):
        for S in itertools.combinations(range(N), k):
            Pi = sub_sym(N, d, S)
            Ck = math.comb(N, k)
            X = pt(P, N, d, S) - Pi / Ck
            m1 = np.linalg.eigvalsh((X + X.T) / 2).min()
            wp, Vp = np.linalg.eigh(Pi)
            B = Vp[:, wp > 0.5]
            lam = np.linalg.eigvalsh(B.T @ pt(P / D, N, d, S) @ B).min()
            if m1 < -1e-8 or abs(lam - 1 / (D * Ck)) > 1e-8:
                ok = False; print("FAIL(1)", N, d, S, m1, lam, 1 / (D * Ck))
    # (2), (3)
    for trial in range(30):
        z = rng.normal(size=D) + 1j * rng.normal(size=D)
        psi = symbasis @ z; psi /= np.linalg.norm(psi)
        Rpsi = np.outer(psi, psi.conj())
        for k in range(1, N):
            for S in itertools.combinations(range(N), k):
                m2 = np.linalg.eigvalsh(pt(Rpsi, N, d, S)).min()
                rho = pmin / D * P + (1 - pmin) * Rpsi
                m3 = np.linalg.eigvalsh(pt(rho, N, d, S)).min()
                if m2 < -0.5 - 1e-9 or m3 < -1e-9:
                    ok = False; print("FAIL(2/3)", N, d, S, m2, m3)
    # GHZ-type extremal state: (|0..0> - |1..1>)/sqrt 2
    g = np.zeros(d ** N, dtype=complex); g[0] = 1; g[np.ravel_multi_index((1,) * N, (d,) * N)] = -1; g /= math.sqrt(2)
    Rg = np.outer(g, g.conj())
    S = tuple(range(N // 2))
    at = np.linalg.eigvalsh(pt(pmin / D * P + (1 - pmin) * Rg, N, d, S)).min()
    p2 = pmin - 1e-3
    below = np.linalg.eigvalsh(pt(p2 / D * P + (1 - p2) * Rg, N, d, S)).min()
    if at < -1e-9 or below >= 0:
        ok = False; print("FAIL(GHZ)", N, d, at, below)
    print(f"N={N} d={d} D={D} pmin={pmin:.6f} ghz_at_pmin={at:.2e} ghz_below={below:.2e}")
print("RESULT", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
