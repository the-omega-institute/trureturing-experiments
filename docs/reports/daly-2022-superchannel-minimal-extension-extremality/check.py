"""Daly, "Axiomatic Approach to Quantum Superchannels", arXiv:2210.00370v1, Remark 4.1: are the extensions of a QSC
that give minimal dimension e = rank(Tr_{r1} Tr_{r2} C_S~) extreme points of the extension set?

Search-seat refutation with d1 = r1 = d2 = r2 = 2.  Channel Choi matrices X = (phi(E_ij)) live in M_A (x) M_B (A = d1,
B = r1, index (a,b) a-major); outputs in M_D (x) M_C (D = d2, C = r2).
  F_s(X) = X/2 + (I_D (x) Tr_A X)/8 + Tr(X) I/16 + s Tr[(Z_A (x) I_B) X] (Z_D (x) I_C)/32,  -1 <= s <= 1.
Checks (exact, sympy, unless marked float):
 1. S(2,2) = {X : Tr_B X = t I_A} has dimension 13, and F_s - F_0 vanishes on a basis of S(2,2).
 2. The superchannel Choi matrix J_s = (F_s(E_IJ))_{IJ} is positive definite for s = -1, 0, 1 (all eigenvalues).
 3. F_s maps Tr_B X = I_A to Tr_C F_s(X) = I_D (TP maps to TP maps) for symbolic X.
 4. R_s = Tr_{r1} Tr_{r2} J_s = P + I/2 + (s/8) Z (x) Z with eigenvalues 5/2 + s/8, 1/2 + s/8, 1/2 - s/8 (twice).
 5. For b != b', the Choi block of J_0 between input rows (a,b) and columns (a',b') is K (x) E^C_{bb'} on (d,c),
    with K = P/2 + I/8 positive definite; these input matrix units lie in S(2,2).
 6. F_1 != F_{-1} and F_0 = (F_1 + F_{-1})/2.
 7. (float) 300 random directions G - F_0 vanishing on S(2,2), each scaled to the boundary of the CP cone: the reduced
    matrix R of the boundary extension stays nonsingular (Lemma 2 of the seat: every extension has e = 4).
"""
import random
import sympy as sp
import numpy as np

R2 = range(2)
Z = sp.diag(1, -1); I2 = sp.eye(2)
def idx(a, b): return 2 * a + b
def E(n, i, j):
    M = sp.zeros(n, n); M[i, j] = 1; return M
def trA(X): return sp.Matrix(2, 2, lambda b, bb: sum(X[idx(a, b), idx(a, bb)] for a in R2))
def trB(X): return sp.Matrix(2, 2, lambda a, aa: sum(X[idx(a, b), idx(aa, b)] for b in R2))
def trC(Y): return sp.Matrix(2, 2, lambda d, dd: sum(Y[idx(d, c), idx(dd, c)] for c in R2))
ZI = sp.kronecker_product(Z, I2)
def F(X, s):
    return X / 2 + sp.kronecker_product(I2, trA(X)) / 8 + X.trace() * sp.eye(4) / 16 + s * (ZI * X).trace() * ZI / 32
def choi_of(Fn):
    M = sp.zeros(16, 16)
    for I in range(4):
        for J in range(4):
            Y = Fn(E(4, I, J))
            for K in range(4):
                for L in range(4):
                    M[4 * I + K, 4 * J + L] = Y[K, L]
    return M
def red(M):  # rows ((a,b),(d,c)) -> (a,d)
    return sp.Matrix(4, 4, lambda p, q: sum(M[4 * idx(p // 2, b) + idx(p % 2, c), 4 * idx(q // 2, b) + idx(q % 2, c)]
                                            for b in R2 for c in R2))
ok = True
# 1. S(2,2) basis
basis = []
for I in range(4):
    for J in range(4):
        a, b = divmod(I, 2); aa, bb = divmod(J, 2)
        X = E(4, I, J)
        if b != bb:
            basis.append(X)
for a in R2:
    for aa in R2:
        basis.append(E(4, idx(a, 0), idx(aa, 0)) - E(4, idx(a, 1), idx(aa, 1)))
basis.append(sp.eye(4))
inS = all(trB(X)[0, 1] == 0 and trB(X)[1, 0] == 0 and trB(X)[0, 0] == trB(X)[1, 1] for X in basis)
dimS = sp.Matrix([list(X) for X in basis]).rank()
vanish = all((F(X, 1) - F(X, 0)).is_zero_matrix for X in basis)
print("1. basis in S:", inS, "dim S =", dimS, "(expected 13); F_s - F_0 vanishes on S:", vanish)
ok &= inS and dimS == 13 and vanish
# 2. positivity
for sv in (-1, 0, 1):
    ev = choi_of(lambda X: F(X, sv)).eigenvals()
    mn = min(ev)
    print(f"2. s={sv}: J_s eigenvalues {dict(sorted(ev.items()))}, minimum {mn}")
    ok &= mn > 0
# 3. TP preservation
x = sp.symbols('x0:16')
X = sp.Matrix(4, 4, x)
sub = {X[idx(0, 1), idx(0, 1)]: 1 - X[idx(0, 0), idx(0, 0)], X[idx(1, 1), idx(1, 1)]: 1 - X[idx(1, 0), idx(1, 0)],
       X[idx(0, 1), idx(1, 1)]: -X[idx(0, 0), idx(1, 0)], X[idx(1, 1), idx(0, 1)]: -X[idx(1, 0), idx(0, 0)]}
s = sp.symbols('s')
tp = sp.simplify(trC(F(X, s)).subs(sub)) == sp.eye(2)
print("3. Tr_B X = I  =>  Tr_C F_s(X) = I:", tp)
ok &= tp
# 4. reduced matrix
Rs = sp.simplify(red(choi_of(lambda Y: F(Y, s))))
Om = sp.Matrix([1, 0, 0, 1]); P = Om * Om.T
expectR = P + sp.eye(4) / 2 + s / 8 * sp.kronecker_product(Z, Z)
lam = sp.symbols('lam')
cp = sp.factor(Rs.charpoly(lam).as_expr())
print("4. R_s == P + I/2 + (s/8) Z(x)Z:", sp.simplify(Rs - expectR).is_zero_matrix, "; charpoly", cp)
ok &= sp.simplify(Rs - expectR).is_zero_matrix
ok &= sp.expand(Rs.charpoly(lam).as_expr() - (lam - sp.Rational(5, 2) - s / 8) * (lam - sp.Rational(1, 2) - s / 8)
                * (lam - sp.Rational(1, 2) + s / 8) ** 2) == 0
# 5. coherence blocks
J0 = choi_of(lambda Y: F(Y, 0))
K = P / 2 + sp.eye(4) / 8
blk = True
for b in R2:
    for bb in R2:
        if b == bb: continue
        for a in R2:
            for aa in R2:
                for d in R2:
                    for dd in R2:
                        for c in R2:
                            for cc in R2:
                                v = J0[4 * idx(a, b) + idx(d, c), 4 * idx(aa, bb) + idx(dd, cc)]
                                w = K[idx(a, d), idx(aa, dd)] * (1 if (c, cc) == (b, bb) else 0)
                                blk &= v == w
print("5. off-diagonal blocks equal K (x) E_{bb'}:", blk, "; K eigenvalues", K.eigenvals())
ok &= blk and min(K.eigenvals()) > 0
# 6. distinct and midpoint
d = choi_of(lambda Y: F(Y, 1)) - choi_of(lambda Y: F(Y, -1))
mid = (choi_of(lambda Y: F(Y, 1)) + choi_of(lambda Y: F(Y, -1))) / 2 - J0
print("6. F_1 != F_-1:", not d.is_zero_matrix, "; F_0 = midpoint:", mid.is_zero_matrix, "; delta(Z(x)I) =",
      list((F(ZI, 1) - F(ZI, 0)).diagonal()))
ok &= (not d.is_zero_matrix) and mid.is_zero_matrix
# 7. random boundary extensions (float)
Bmat = np.array([[complex(v) for v in Xb] for Xb in basis])  # 13 x 16, rows = vec(S basis)
_, _, Vh = np.linalg.svd(Bmat.conj())
null = Vh[13:].conj()  # functionals vanishing on S: f(X) = sum f_ij X_ij  (rows span S^perp under bilinear pairing)
# check: f(X) = null_row . vec(X) = 0 for X in S
assert np.allclose(Bmat @ null.T, 0)
J0f = np.array(J0.tolist(), dtype=complex)
rng = np.random.default_rng(221000370)
def red_np(M):
    R = np.zeros((4, 4), dtype=complex)
    for p in range(4):
        for q in range(4):
            R[p, q] = sum(M[4 * idx(p // 2, b) + idx(p % 2, c), 4 * idx(q // 2, b) + idx(q % 2, c)] for b in R2 for c in R2)
    return R
mins = []
for t in range(300):
    # G = F_0 + D with D(X) = sum_k f_k(X) Y_k, Y_k random Hermitian-preserving outputs; Choi of D on E_IJ
    Ys = [rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)) for _ in range(3)]
    JD = np.zeros((16, 16), dtype=complex)
    for I in range(4):
        for Jx in range(4):
            out = sum(null[k][4 * I + Jx] * Ys[k] for k in range(3))
            JD[4 * I:4 * I + 4, 4 * Jx:4 * Jx + 4] = out
    JD = (JD + JD.conj().T) / 2  # Hermitian part keeps the map Hermiticity-preserving and still vanishes on S*-closed S
    lo, hi = 0.0, 1.0
    while np.linalg.eigvalsh(J0f + hi * JD).min() > 0: hi *= 2
    for _ in range(60):
        midt = (lo + hi) / 2
        if np.linalg.eigvalsh(J0f + midt * JD).min() > 0: lo = midt
        else: hi = midt
    Jb = J0f + lo * JD
    mins.append(np.linalg.eigvalsh(red_np(Jb)).min())
print(f"7. (float) 300 boundary extensions: min eigenvalue of R over samples = {min(mins):.6f} (> 0 means e = 4)")
ok &= min(mins) > 1e-6
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
