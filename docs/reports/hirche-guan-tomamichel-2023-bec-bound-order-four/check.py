"""Numerical check of Hirche-Guan-Tomamichel Conjecture V.6 (arXiv:2305.02589v1, eq. V.16) at alpha = 4.

H4(omega) = -(1/3) log tr[((I x r^{-3/8}) omega (I x r^{-3/8}))^4] for a cq state omega with quantum marginal r
(negative powers on the support). K = exp(-3 H4). Inputs rho_i = (|0><0| x s_i0 + |1><1| x s_i1)/2, XOR output
blocks (s10 x s20 + s11 x s21)/4, (s10 x s21 + s11 x s20)/4. Claim: K_out >= (8 K1 K2 - K1 - K2 + 1)/7, and the
exact defect K_out - (8K1K2-K1-K2+1)/7 = [8(a1-b1)(a2-b2) + 4(a1-c1)(a2-c2) + 2(b1-c1)(b2-c2)]/56.
"""
import numpy as np
rng = np.random.default_rng(4)
def rand_density(d, rank=None):
    r = rank or d
    G = rng.normal(size=(d, r)) + 1j * rng.normal(size=(d, r))
    M = G @ G.conj().T
    return M / np.trace(M).real
def mpow(M, p):
    w, V = np.linalg.eigh((M + M.conj().T) / 2)
    out = np.zeros_like(w)
    mask = w > 1e-12
    out[mask] = w[mask] ** p
    return (V * out) @ V.conj().T
def K_of(blocks):
    r = sum(blocks)
    R = mpow(r, -3 / 8)
    return sum(np.trace(np.linalg.matrix_power(R @ B @ R, 4)).real for B in blocks)
def moments(s0, s1):
    r = (s0 + s1) / 2; R = mpow(r, -3 / 8)
    A = R @ s0 @ R / 2; B = R @ s1 @ R / 2; S = A + B; D = A - B
    tr = lambda X: np.trace(X).real
    return tr(S @ S @ D @ D), tr(S @ D @ S @ D), tr(D @ D @ D @ D), tr(np.linalg.matrix_power(S, 4))
bad = 0; worst = np.inf; n = 0
for trial in range(4000):
    d1, d2 = rng.integers(1, 5, size=2)
    rk = lambda d: int(rng.integers(1, d + 1)) if rng.random() < 0.3 else None
    s10, s11 = rand_density(d1, rk(d1)), rand_density(d1, rk(d1))
    s20, s21 = rand_density(d2, rk(d2)), rand_density(d2, rk(d2))
    K1 = K_of([s10 / 2, s11 / 2]); K2 = K_of([s20 / 2, s21 / 2])
    O0 = (np.kron(s10, s20) + np.kron(s11, s21)) / 4; O1 = (np.kron(s10, s21) + np.kron(s11, s20)) / 4
    Ko = K_of([O0, O1])
    bound = (8 * K1 * K2 - K1 - K2 + 1) / 7
    a1, b1, c1, t1 = moments(s10, s11); a2, b2, c2, t2 = moments(s20, s21)
    defect = (8 * (a1 - b1) * (a2 - b2) + 4 * (a1 - c1) * (a2 - c2) + 2 * (b1 - c1) * (b2 - c2)) / 56
    n += 1; worst = min(worst, Ko - bound)
    if Ko - bound < -1e-10: bad += 1; print("VIOLATION", trial, Ko - bound)
    if abs((Ko - bound) - defect) > 1e-9: bad += 1; print("DEFECT MISMATCH", trial, Ko - bound, defect)
    if not (a1 >= b1 - 1e-10 and b1 >= c1 - 1e-10 and c1 >= -1e-12): bad += 1; print("ORDER FAIL", a1, b1, c1)
    if abs(t1 - 1) > 1e-9 and np.linalg.matrix_rank(s10 + s11) == d1: bad += 1; print("TRACE S^4 != 1", t1)
# equality family: erasure sector C plus orthogonal signals P, Q (one input), arbitrary other input
for trial in range(200):
    c = rng.random(); C = c * rand_density(1)
    s10 = np.diag([c, 1 - c, 0]).astype(complex); s11 = np.diag([c, 0, 1 - c]).astype(complex)
    d2 = int(rng.integers(1, 4)); s20, s21 = rand_density(d2), rand_density(d2)
    K1 = K_of([s10 / 2, s11 / 2]); K2 = K_of([s20 / 2, s21 / 2])
    Ko = K_of([(np.kron(s10, s20) + np.kron(s11, s21)) / 4, (np.kron(s10, s21) + np.kron(s11, s20)) / 4])
    if abs(Ko - (8 * K1 * K2 - K1 - K2 + 1) / 7) > 1e-9: bad += 1; print("EQUALITY FAIL", Ko)
print("random trials", n, "min K_out - bound", worst)
print("bad=", bad)
raise SystemExit(1 if bad else 0)
