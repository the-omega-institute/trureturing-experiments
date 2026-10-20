"""Samanta-Mondal-Sen, arXiv:2505.18539v1, Appendix B.3: E_N(U) = D_N(U) for diagonal unitaries U.

The proof uses coordinatewise complex conjugation C in the computational basis: C maps fully separable and biseparable
pure states onto themselves, preserves |<phi|psi>|^2, and U^dagger C = C U for diagonal U. This script checks the
pointwise identity G(U p) = G(U^dagger C p) for random diagonal phases and random product inputs, N = 3, 4, 5, with G the
generalized geometric measure computed from the largest Schmidt coefficient over all nontrivial bipartitions
(the source's Appendix A.1 formula).
"""
import itertools
import sys
import numpy as np

rng = np.random.default_rng(7)


def ggm(v, N):
    psi = v.reshape([2] * N)
    best = 0.0
    for r in range(1, N // 2 + 1):
        for A in itertools.combinations(range(N), r):
            B = [i for i in range(N) if i not in A]
            M = np.transpose(psi, list(A) + B).reshape(2 ** len(A), -1)
            best = max(best, np.linalg.svd(M, compute_uv=False)[0] ** 2)
    return 1 - best


def product(ps):
    v = np.array([1.0 + 0j])
    for p in ps:
        v = np.kron(v, p)
    return v


maxdiff = 0.0
cases = 0
for N in (3, 4, 5):
    for _ in range(200):
        U = np.exp(1j * rng.uniform(0, 2 * np.pi, 2 ** N))
        ps = []
        for _i in range(N):
            a = rng.normal(size=2) + 1j * rng.normal(size=2)
            ps.append(a / np.linalg.norm(a))
        p = product(ps)
        maxdiff = max(maxdiff, abs(ggm(U * p, N) - ggm(np.conj(U) * np.conj(p), N)))
        cases += 1
print(f"cases: {cases}; max |G(U p) - G(U^dagger C p)| = {maxdiff}")
ok = maxdiff < 1e-10
print("RESULT PASS" if ok else "RESULT FAIL")
sys.exit(0 if ok else 1)
