"""Cheng-Chen, arXiv:2605.19530v2, Conjecture 24: for odd n, the range of I_rho = Tr(rho^T eps_n rho eps_n) over
n-qubit states rho is [-1/2, 0], with eps = [[0, 1], [-1, 0]] and eps_n = eps^{tensor n}.

Checks (numpy, double precision, plus exact endpoint evaluations with Fractions):
  - for odd n in {1, 3, 5}: random states of every rank have I real and in [-1/2, 0] (up to 1e-12);
  - eps_n^T = -eps_n and eps_n^2 = -Id for odd n; for even n the signs are + (for contrast);
  - I(pure state) = 0 and I((|v><v| + |Theta v><Theta v|)/2) = -1/2, Theta v = eps_n conj(v);
  - the intermediate value: I along the segment between those two states covers [-1/2, 0];
  - for even n the range is different (an example with I > 0).
Usage: python3 check.py SAMPLES SEED
"""
import sys
import numpy as np
from fractions import Fraction

def eps_n(n):
    e = np.array([[0, 1], [-1, 0]], dtype=complex); m = np.array([[1]], dtype=complex)
    for _ in range(n): m = np.kron(m, e)
    return m

def I(rho, E): return np.trace(rho.T @ E @ rho @ E)

def rand_state(d, rank, rng):
    G = rng.normal(size=(d, rank)) + 1j * rng.normal(size=(d, rank))
    r = G @ G.conj().T; return r / np.trace(r)

def main():
    S, seed = int(sys.argv[1]), int(sys.argv[2]); rng = np.random.default_rng(seed); ok = True
    for n in [1, 3, 5]:
        E = eps_n(n); d = 2 ** n
        sym = np.allclose(E.T, -E) and np.allclose(E @ E, -np.eye(d))
        vals = []
        for _ in range(S):
            rho = rand_state(d, int(rng.integers(1, d + 1)), rng); v = I(rho, E); vals.append(v)
        re = np.array([x.real for x in vals]); im = max(abs(x.imag) for x in vals)
        inrange = re.min() >= -0.5 - 1e-12 and re.max() <= 1e-12 and im < 1e-12
        v = rng.normal(size=d) + 1j * rng.normal(size=d); v /= np.linalg.norm(v); tv = E @ v.conj()
        P = np.outer(v, v.conj()); Q = np.outer(tv, tv.conj())
        e0 = I(P, E); e1 = I((P + Q) / 2, E)
        seg = [I((1 - t) * P + t * (P + Q) / 2, E).real for t in np.linspace(0, 1, 101)]
        cover = abs(seg[0]) < 1e-12 and abs(seg[-1] + 0.5) < 1e-12 and all(seg[i + 1] <= seg[i] + 1e-12 for i in range(100))
        print(f"n={n}: eps^T=-eps, eps^2=-I: {sym}; {S} random states in [-1/2,0] (min {re.min():.6f}, max {re.max():.2e}, |Im| < {im:.1e}): {inrange}; pure {e0.real:.2e}; mixed pair {e1.real:.12f}; segment monotone from 0 to -1/2: {cover}")
        ok &= sym and inrange and abs(e0) < 1e-12 and abs(e1 + 0.5) < 1e-12 and cover
    # exact endpoints at n = 1 with Fractions: rho = |0><0| gives 0, rho = I/2 gives -1/2
    def mul(A, B): return [[sum(A[i][k] * B[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
    e = [[0, 1], [-1, 0]]; T = lambda A: [[A[j][i] for j in range(2)] for i in range(2)]
    tr = lambda A: A[0][0] + A[1][1]
    r0 = [[Fraction(1), 0], [0, 0]]; rh = [[Fraction(1, 2), 0], [0, Fraction(1, 2)]]
    x0 = tr(mul(mul(mul(T(r0), e), r0), e)); xh = tr(mul(mul(mul(T(rh), e), rh), e))
    print(f"exact n=1: I(|0><0|) = {x0}, I(Id/2) = {xh}"); ok &= x0 == 0 and xh == Fraction(-1, 2)
    E2 = eps_n(2); bell = np.zeros(4, dtype=complex); bell[1] = bell[2] = 1 / np.sqrt(2)
    x2 = I(np.outer(bell, bell.conj()), E2).real
    print(f"contrast n=2: eps^T = +eps: {np.allclose(E2.T, E2)}; I at a Bell state = {x2:.6f} (outside [-1/2, 0])"); ok &= x2 > 0
    print("RESULT", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)

main()
