"""Puchała–Jenčová–Sedlák–Ziman arXiv:1504.00477v1, Section V: the conjecture b(E (x) F) = b(E) b(F) for channels fails.
Phi_t(X) = sum_ij P_t[i][j] X_jj |i><i|, P_t = [[1/5,9/20,1/4],[1/2,t,3/10],[3/10,11/20-t,9/20]], 0 < t <= 1/10.
1. Exact (Fractions): column sums; the dual certificate r, c with r_i + c_j >= 1/P_ij; max over the 6 permutations of
   sum 1/P_{i,pi(i)} equals m(t) = 1/t + 22/3; for the tensor product, the best of all 9! permutations of 1/(P_t (x) P_s)
   is at least m(t) m(s) + 229/81 (the swap sigma), checked for t, s in {1/10, 1/20, 1/50}.
2. Definition check (floats): b(Phi) = inf over channels Psi of t_Phi(Psi) = 1/lambda_max(J^-1/2 J_Psi J^-1/2) (capped at 1);
   random channels and random unitary channels never go below 1/m(t); the permutation unitary pi attains it; for the
   tensor product the unitary channel of sigma gives t = 1/(m^2 + 229/81) < b(Phi)^2.
"""
import itertools, math, sys
from fractions import Fraction as F
import numpy as np

rng = np.random.default_rng(20261017)

def P(t): return [[F(1, 5), F(9, 20), F(1, 4)], [F(1, 2), t, F(3, 10)], [F(3, 10), F(11, 20) - t, F(9, 20)]]
def m(t): return 1 / t + F(22, 3)

def assign_max(A):
    n = len(A); return max(sum(A[i][p[i]] for i in range(n)) for p in itertools.permutations(range(n)))

def kron(A, B):
    n, k = len(A), len(B)
    return [[A[i // k][j // k] * B[i % k][j % k] for j in range(n * k)] for i in range(n * k)]

def choi_diag(Pm):  # normalized Choi of the measure-and-prepare channel: diag(P_ij)/d, index (i out, j in)
    d = len(Pm); return np.diag([float(Pm[i][j]) / d for i in range(d) for j in range(d)])

def choi_of_kraus(Ks, d):
    J = np.zeros((d * d, d * d), complex)
    for j in range(d):
        for jj in range(d):
            Ejj = np.zeros((d, d)); Ejj[j, jj] = 1
            out = sum(K @ Ejj @ K.conj().T for K in Ks) / d
            for i in range(d):
                for ii in range(d):
                    J[i * d + j, ii * d + jj] = out[i, ii]
    return J

def t_val(JPhi, JPsi):
    w, V = np.linalg.eigh(JPhi); R = V @ np.diag(w ** -0.5) @ V.conj().T
    lam = np.linalg.eigvalsh(R @ JPsi @ R).max()
    return min(1.0, 1 / lam)

def main():
    fails = []
    for t in (F(1, 10), F(1, 20), F(1, 50)):
        Pt = P(t)
        if any(sum(Pt[i][j] for i in range(3)) != 1 for j in range(3)) or min(min(r) for r in Pt) <= 0: fails.append(("stochastic", t))
        A = [[1 / x for x in row] for row in Pt]
        r = [F(3), F(7, 3), F(4, 3)]; c = [F(2), 1 / t - F(7, 3), F(1)]
        if any(r[i] + c[j] < A[i][j] for i in range(3) for j in range(3)) or sum(r) + sum(c) != m(t): fails.append(("dual", t))
        if assign_max(A) != m(t): fails.append(("onecopy", t))
    print("one-copy: column-stochastic, dual certificate and exact optimum m(t) = 1/t + 22/3 for t in {1/10, 1/20, 1/50}:", not fails)
    for t, s in ((F(1, 20), F(1, 20)), (F(1, 10), F(1, 20)), (F(1, 50), F(1, 10))):
        At = [[1 / x for x in row] for row in P(t)]; As = [[1 / x for x in row] for row in P(s)]
        K = kron(At, As)
        best = assign_max(K)
        sigma = [0, 7, 6, 5, 4, 3, 2, 1, 8]
        sig = sum(K[i][sigma[i]] for i in range(9))
        ok = sig == m(t) * m(s) + F(229, 81) and best >= sig
        print(f"t={t}, s={s}: m(t)m(s) = {m(t)*m(s)}, sigma objective = {sig} (= +229/81: {sig == m(t)*m(s) + F(229,81)}), best over 9! = {best}; "
              f"b(Phi_t)b(Phi_s) = {1/(m(t)*m(s))} > b(tensor) <= {1/best}")
        if not ok: fails.append(("tensor", t, s))
    # definition check
    t = F(1, 20); Pt = P(t); J = choi_diag(Pt); target = 1 / float(m(t))
    worst = 1.0
    for trial in range(4000):
        if trial % 2:
            Z = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); U, _ = np.linalg.qr(Z); Ks = [U]
        else:
            k = int(rng.integers(1, 6)); G = rng.normal(size=(3 * k, 3)) + 1j * rng.normal(size=(3 * k, 3)); Q, _ = np.linalg.qr(G)
            Ks = [Q[3 * a:3 * a + 3, :] for a in range(k)]
        worst = min(worst, t_val(J, choi_of_kraus(Ks, 3)))
    perm = np.zeros((3, 3)); perm[2, 0] = perm[1, 1] = perm[0, 2] = 1  # pi: column j -> row pi^-1 ; maps |j> to |2-j>
    tperm = t_val(J, choi_of_kraus([perm], 3))
    print(f"definition check (t = 1/20): min sampled t_Phi(Psi) = {worst:.6f} >= 1/m = {target:.6f}; permutation unitary gives {tperm:.6f}")
    if worst < target - 1e-9 or abs(tperm - target) > 1e-9: fails.append("definition")
    J2 = np.kron(J, J) * 1.0
    # reorder Kronecker of Choi matrices to the joint (out1 out2, in1 in2) ordering via the measure-and-prepare form
    Pk = kron(Pt, Pt); J2 = choi_diag(Pk)
    S = np.zeros((9, 9)); sigma = [0, 7, 6, 5, 4, 3, 2, 1, 8]
    for i in range(9): S[i, sigma[i]] = 1    # row i -> column sigma(i): the unitary has entry (i, sigma(i))
    tsig = t_val(J2, choi_of_kraus([S], 9))
    bound = 1 / float(m(t) ** 2 + F(229, 81))
    print(f"tensor definition check: unitary sigma gives t = {tsig:.8f} = 1/(m^2 + 229/81) = {bound:.8f} < b(Phi)^2 = {target**2:.8f}")
    if abs(tsig - bound) > 1e-9 or not tsig < target ** 2: fails.append("tensor-definition")
    print("RESULT", "PASS" if not fails else f"FAIL {fails}")
    sys.exit(0 if not fails else 1)

if __name__ == "__main__":
    main()
