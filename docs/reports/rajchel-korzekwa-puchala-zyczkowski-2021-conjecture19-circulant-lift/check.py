#!/usr/bin/env python3
"""Circulant unistochastic matrices without a circulant unitary lift (Conjecture 19 refuted).

Source: Rajchel-Mieldzioc, Korzekwa, Puchala, Zyczkowski, arXiv:2101.11288; J. Math. Phys. 63, 012202 (2022),
Conjecture 19: "The sets C_d cap U_d and C_d^2 coincide for all d, i.e., for every circulant unistochastic matrix
it is possible to find a corresponding circulant unitary matrix."

For N = 2m (m >= 5), S = {0, 1, N-1, m-1, m+1}, T the m-cycle shift, R = [[1,1],[-1,-1]],
Q = R (x) T - R^T (x) T^T and U = (3 I + 2 Q)/5.  Exact integer arithmetic (scaled by 5) checks:
  * U U^T = I, so B = (U_ij^2) is orthostochastic (hence unistochastic);
  * B is circulant with first row b_0 = 9/25, b_s = 4/25 on {1, N-1, m-1, m+1}, 0 elsewhere;
  * the support overlaps used by the phase argument: for h = 1, m-1, m-2 the pairs (s, s-h) with both in S
    are exactly {(0,N-1),(1,0)}, {(0,m+1),(m-1,0)}, {(m-1,1),(N-1,m+1)}.
Phase argument (proof, not computed): a circulant unitary V = circ(z) with |z_s|^2 = b_s has z_s != 0 exactly on S;
orthogonality at h = 1 and h = m-1 gives z_1 = -u conj(z_{N-1}) and z_{m-1} = -u conj(z_{m+1}) with u = z_0/conj(z_0);
at h = m-2 it gives 0 = z_{N-1} conj(z_{m+1}) + z_{m-1} conj(z_1) = 2 z_{N-1} conj(z_{m+1}), a contradiction.
Prints one line per N and ALL_OK.
"""
import sys
def fail(m): print("FAIL", m); sys.exit(1)
for m in range(5, 16):
    N = 2 * m
    T = [[1 if c == (r + 1) % m else 0 for c in range(m)] for r in range(m)]
    Rm = [[1, 1], [-1, -1]]
    def kron(A, B):
        return [[A[i // len(B)][j // len(B[0])] * B[i % len(B)][j % len(B[0])] for j in range(len(A[0]) * len(B[0]))] for i in range(len(A) * len(B))]
    def tr(A): return [list(x) for x in zip(*A)]
    K1 = kron(Rm, T); K2 = kron(tr(Rm), tr(T))
    Q = [[K1[i][j] - K2[i][j] for j in range(N)] for i in range(N)]
    U5 = [[(3 if i == j else 0) + 2 * Q[i][j] for j in range(N)] for i in range(N)]   # 5 U
    for i in range(N):
        for j in range(N):
            if sum(U5[i][k] * U5[j][k] for k in range(N)) != (25 if i == j else 0): fail(f"orthogonal N={N}")
    S = {0, 1, N - 1, m - 1, m + 1}
    row0 = [U5[0][j] ** 2 for j in range(N)]
    want = [9 if s == 0 else (4 if s in S else 0) for s in range(N)]
    if row0 != want: fail(f"row N={N}")
    for i in range(N):
        if [U5[i][j] ** 2 for j in range(N)] != [row0[(j - i) % N] for j in range(N)]: fail(f"circulant N={N}")
    pairs = {h: sorted((s, (s - h) % N) for s in S if (s - h) % N in S) for h in (1, m - 1, m - 2)}
    if pairs[1] != sorted([(0, N - 1), (1, 0)]): fail(f"h=1 N={N}")
    if pairs[m - 1] != sorted([(0, m + 1), (m - 1, 0)]): fail(f"h=m-1 N={N}")
    if pairs[m - 2] != sorted([(m - 1, 1), (N - 1, m + 1)]): fail(f"h=m-2 N={N}")
    print(f"N={N}: U orthogonal, B circulant with first row (x25) {row0}, overlaps {pairs}")
print("ALL_OK")
