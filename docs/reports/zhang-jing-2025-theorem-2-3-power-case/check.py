"""Zhang-Jing, ELA 41 (2025), the open case of Theorem 2.3: does
    [[|A|^{2a}, A*], [A, |A|^{2(1-a)}]] >= 0  iff  A is normal
hold for 1/2 <= a <= 1?  Expected: yes for 1/2 <= a < 1 (every A, singular included), no at a = 1."""
import sys
import numpy as np
from fractions import Fraction as Fr

ok = True
def check(name, cond):
    global ok
    print(("PASS " if cond else "FAIL ") + name); ok = ok and bool(cond)

def psd_pow(M, e):
    w, V = np.linalg.eigh((M + M.conj().T) / 2); w = np.clip(w, 0, None)
    p = np.where(w > 1e-13, w ** e, 0.0) if e > 0 else np.ones_like(w)   # x^0 = 1
    return (V * p) @ V.conj().T
def block(A, a):
    D = psd_pow(A.conj().T @ A, 0.5)
    return np.block([[psd_pow(D @ D, a), A.conj().T], [A, psd_pow(D @ D, 1 - a)]])
def mineig(M): return np.linalg.eigvalsh((M + M.conj().T) / 2).min()
def normal(A): return np.allclose(A.conj().T @ A, A @ A.conj().T, atol=1e-9)

rng = np.random.default_rng(20261018)
def rand(n, rank=None, normal_=False):
    if normal_:
        Q, _ = np.linalg.qr(rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n)))
        lam = rng.standard_normal(n) + 1j * rng.standard_normal(n)
        if rank is not None: lam[rank:] = 0
        return Q @ np.diag(lam) @ Q.conj().T
    X = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    if rank is not None:
        U, s, Vh = np.linalg.svd(X); s[rank:] = 0; X = (U * s) @ Vh
    return X

grid = [0.5, 0.55, 0.6, 0.75, 0.9, 0.99]
bad_nonnormal = 0; bad_normal = 0; trials = 0
for n in (2, 3, 4):
    for rank in (None, n - 1, 1):
        for t in range(25):
            for a in grid:
                A = rand(n, rank); An = rand(n, rank, normal_=True); trials += 1
                if not normal(A) and mineig(block(A, a)) > -1e-9: bad_nonnormal += 1
                if mineig(block(An, a)) < -1e-8: bad_normal += 1
check(f"1/2 <= a < 1: every random nonnormal A (singular included) gives a non-PSD block ({trials} pairs, failures {bad_nonnormal})", bad_nonnormal == 0)
check(f"1/2 <= a < 1: every random normal A (singular included) gives a PSD block (failures {bad_normal})", bad_normal == 0)

A1 = np.array([[0, 1], [2, 0]], dtype=complex)
check(f"a = 1: A = [[0, 1], [2, 0]] is nonnormal and the block is PSD (min eigenvalue {mineig(block(A1, 1.0)):.2e})",
      not normal(A1) and mineig(block(A1, 1.0)) > -1e-12)
# exact a = 1 certificate: [[D^2, A*], [A, I]] = [A*; I] [A, I] restricted, i.e. the block equals Q* Q with Q = [A, I]
Q = np.hstack([A1, np.eye(2)])
check("a = 1: the block equals Q*Q with Q = [A, I] (so it is PSD for every A)", np.allclose(block(A1, 1.0), Q.conj().T @ Q))
# exact spot check of the mechanism at a = 3/4 for a singular nonnormal A
As = np.array([[0, 1], [0, 0]], dtype=complex)
check(f"a = 3/4, singular nonnormal A = [[0, 1], [0, 0]]: block not PSD (min eigenvalue {mineig(block(As, 0.75)):.4f})", mineig(block(As, 0.75)) < -1e-6)
print("RESULT", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
