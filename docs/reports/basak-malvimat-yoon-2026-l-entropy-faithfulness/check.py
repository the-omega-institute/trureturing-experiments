"""Check of a counterexample to the L-entropy faithfulness conjecture of Basak-Malvimat-Yoon,
arXiv:2602.00617v1 (PRL 136, 080202 (2026)), p. 2: "We conjecture that it remains nonzero
for any GME state".

l_{ij} = min{2 S(A_i), 2 S(A_j)} - S_R(A_i : A_j), S_R = S(A A*) in |sqrt(rho_AB)>, and
l_all = (prod_{i<j} l_{ij})^{2/(n(n-1))}.  Entropies in bits.
Witness: three parties A = C^4 = A_L (x) A_R, B = C^2, C = C^2,
|Psi> = |Phi>_{A_L B} (x) |Phi>_{A_R C}.  Family: |phi>_{A_L B} (x) |chi>_{A_R C}.
Prediction registered before the run: the witness is GME (all three marginals of rank > 1),
l_AB = l_AC = 0, l_BC = 2, l_all = 0; every family member with entangled phi, chi is GME with
l_AB = 0; random tripartite states of the same dimensions have l_all > 0.

Usage: python3 check.py [seed]   (numpy)
"""
import sys
import numpy as np


def entropy(rho):
    ev = np.clip(np.linalg.eigvalsh((rho + rho.conj().T) / 2), 0, None)
    ev = ev[ev > 1e-14]
    return float(-(ev * np.log2(ev)).sum())


def ptrace(psi, dims, keep):
    """Reduced density matrix of a pure state on the parties in keep (sorted)."""
    n = len(dims)
    t = psi.reshape(dims)
    rest = [i for i in range(n) if i not in keep]
    t = np.transpose(t, list(keep) + rest)
    dk = int(np.prod([dims[i] for i in keep]))
    m = t.reshape(dk, -1)
    return m @ m.conj().T


def sqrtm_psd(rho):
    ev, U = np.linalg.eigh((rho + rho.conj().T) / 2)
    return (U * np.sqrt(np.clip(ev, 0, None))) @ U.conj().T


def reflected(psi, dims, i, j):
    """S_R(A_i : A_j) = S(A_i A_i*) of the canonical purification |sqrt(rho_{ij})>."""
    di, dj = dims[i], dims[j]
    rho = ptrace(psi, dims, sorted([i, j]))
    s = sqrtm_psd(rho)  # operator on A_i (x) A_j (order i < j)
    # |sqrt rho> = sum s[(a b),(a' b')] |a b a' b'>; A_i A_i* marginal
    t = s.reshape(di, dj, di, dj)  # a, b, a', b'
    m = np.transpose(t, (0, 2, 1, 3)).reshape(di * di, dj * dj)
    return entropy(m @ m.conj().T)


def l_entropies(psi, dims):
    n = len(dims)
    S = [entropy(ptrace(psi, dims, [k])) for k in range(n)]
    out = {}
    for i in range(n):
        for j in range(i + 1, n):
            out[(i, j)] = min(2 * S[i], 2 * S[j]) - reflected(psi, dims, i, j)
    vals = list(out.values())
    lall = float(np.prod(np.clip(vals, 0, None)) ** (2 / (n * (n - 1))))
    return S, out, lall


def gme(psi, dims):
    return all(np.linalg.matrix_rank(ptrace(psi, dims, [k]), tol=1e-10) > 1 for k in range(len(dims)))


def bell(d=2):
    v = np.zeros(d * d); v[[k * d + k for k in range(d)]] = 1 / np.sqrt(d); return v


def two_pair(phi, chi):
    """phi on A_L B (2x2), chi on A_R C (2x2) -> state on A (4) x B (2) x C (2), A = (A_L, A_R)."""
    t = np.einsum('xb,yc->xybc', phi.reshape(2, 2), chi.reshape(2, 2))  # A_L, A_R, B, C
    return t.reshape(4, 2, 2).reshape(-1)


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20261016
    rng = np.random.default_rng(seed)
    dims = [4, 2, 2]
    ok = True
    psi = two_pair(bell(), bell())
    S, l, lall = l_entropies(psi, dims)
    print("witness: GME", gme(psi, dims), " S =", [round(x, 12) for x in S], " l =", {k: round(v, 12) for k, v in l.items()}, " l_all =", lall)
    ok &= gme(psi, dims) and abs(l[(0, 1)]) < 1e-9 and abs(l[(0, 2)]) < 1e-9 and abs(l[(1, 2)] - 2) < 1e-9 and lall < 1e-9
    worst = 0.0
    for t in range(200):
        phi = rng.standard_normal(4) + 1j * rng.standard_normal(4); phi /= np.linalg.norm(phi)
        chi = rng.standard_normal(4) + 1j * rng.standard_normal(4); chi /= np.linalg.norm(chi)
        p = two_pair(phi, chi)
        S, l, lall = l_entropies(p, dims)
        ok &= gme(p, dims)
        worst = max(worst, abs(l[(0, 1)]))
    print("200 random family members: all GME, max |l_AB| =", worst)
    ok &= worst < 1e-9
    mins = []
    for t in range(200):
        p = rng.standard_normal(16) + 1j * rng.standard_normal(16); p /= np.linalg.norm(p)
        mins.append(l_entropies(p, dims)[2])
    print("200 random states of C^4 x C^2 x C^2: min l_all =", min(mins))
    ok &= min(mins) > 1e-6
    print("ALL_OK" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
