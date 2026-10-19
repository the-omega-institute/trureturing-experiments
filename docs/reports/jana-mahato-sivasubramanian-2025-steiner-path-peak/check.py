#!/usr/bin/env python3
"""Jana-Mahato-Sivasubramanian (arXiv:2407.03309v1; Comput. Appl. Math. 44 (2025) 223), Conjecture 29: for the path
P_n (n > 5), with D = D_2(P_n)[B, B] (B = (e_1, f_1, ..., e_{n-2}, f_{n-2}, e_{n-1}), e_i = {i, i+1}, f_j = {j, j+2},
entry = 2-Steiner distance = number of edges of the smallest subpath containing the four vertices) and
det(xI - D) = sum_i a_i x^i, every maximiser l of |a_i| over 0 <= i <= 2n-4 equals n - 1.

Exact integer characteristic polynomials (Faddeev-LeVerrier) for 6 <= n <= 22. Checks the source's facts
|a_{2n-4}| = 3n - 5 and (Conjecture 30, settled) a_{2n-5} = -(n-1)(n-2)(2n^2+6n-15)/6, and the peak bound
l <= floor(7n/5). Reports the maximisers; the conjecture fails exactly when a maximiser differs from n - 1.
Exit 0 and RESULT PASS when the expected pattern (peak n - 1 for n <= 15, peak n for 16 <= n <= 22) holds.
"""
import sys
def Dmat(n):
    B = []
    for i in range(1, n - 1): B.append((i, i + 1)); B.append((i, i + 2))
    B.append((n - 1, n))
    def steiner(S): return max(S) - min(S)   # smallest subpath of P_n containing S
    return [[steiner(B[a] + B[c]) for c in range(len(B))] for a in range(len(B))]
def charpoly(A):
    n = len(A); M = [[0] * n for _ in range(n)]; c = [0] * (n + 1); c[n] = 1
    for k in range(1, n + 1):
        AM = [[sum(A[i][l] * M[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        M = [[AM[i][j] + (c[n - k + 1] if i == j else 0) for j in range(n)] for i in range(n)]
        AMk = [[sum(A[i][l] * M[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        tr = sum(AMk[i][i] for i in range(n)); assert tr % k == 0; c[n - k] = -tr // k
    return c
ok = True
for n in range(6, 23):
    c = charpoly(Dmat(n)); N = 2 * n - 3
    ok &= abs(c[2 * n - 4]) == 3 * n - 5
    ok &= 6 * c[2 * n - 5] == -(n - 1) * (n - 2) * (2 * n * n + 6 * n - 15)
    co = [abs(v) for v in c[:2 * n - 3]]; mx = max(co); ells = [i for i, v in enumerate(co) if v == mx]
    ok &= all(l <= (7 * n) // 5 for l in ells)
    expected = [n - 1] if n <= 15 else [n]
    ok &= ells == expected
    line = f"n={n} size={N} maximisers={ells} conjectured={n - 1}"
    if n == 16: line += f"  |a_15|={co[15]}  |a_16|={co[16]}  |a_17|={co[17]}"
    print(line, flush=True)
print("RESULT", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
