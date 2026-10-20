"""Sobolewski-Ulas, "Solutions of certain meta-Fibonacci recurrences", arXiv:2204.04011v1, Question 7.5: with
h(n) = 1 for n <= 1 and h(n) = h(n - h(n-1)) + h(n-2) for n > 1, r = h mod 2 and u_i the increasing enumeration of
{m >= 0 : r(m) = i}, is u_i k-regular for some k?  Search-seat answer: both are 2-regular.

Checks (exact integers, n < 2^18):
 1. Theorem 3.1: h(2n+1) = n+1 and h(2n+2) = h(2n) + h(n+1).
 2. Theorem 3.9 relations: r(8n+2) = 1 - r(8n), r(8n+6) = r(8n+4), r(8n+4) = 0,1,1,0 periodic, r(16n+8) = 0,0,1,1
    periodic, r(32n) = r(16n), r(32n+16) = 1 - r(16n+8).
 3. Block balance: every aligned block r(16q..16q+15) has exactly eight ones.
 4. Block formula: the block is determined by (r(16q), q mod 4), so u_i(8q+j) = 16q + F_i(r(16q), q mod 4, j) for a fixed
    table F_i (printed).
 5. 2-kernel ranks: the span of {u_i(2^e n + b) : e <= E, b < 2^e} restricted to the first 200 terms has the same rank for
    E = 6, 7, 8 (stabilization; evidence for finite generation).
"""
import numpy as np
N = 1 << 18
h = [1] * N
for n in range(2, N):
    a = n - h[n - 1]
    h[n] = (h[a] if a >= 0 else 1) + h[n - 2]
ok = True
t31 = all(h[2*n+1] == n+1 for n in range(N//2 - 1)) and all(h[2*n+2] == h[2*n] + h[n+1] for n in range(N//2 - 1))
print("1. Theorem 3.1:", t31); ok &= t31
r = [x % 2 for x in h]
t39 = (all(r[8*n+2] == 1 - r[8*n] for n in range(N//8)) and all(r[8*n+6] == r[8*n+4] for n in range(N//8))
       and all(r[8*n+4] == [0,1,1,0][n % 4] for n in range(N//8)) and all(r[16*n+8] == [0,0,1,1][n % 4] for n in range(N//16))
       and all(r[32*n] == r[16*n] for n in range(N//32)) and all(r[32*n+16] == 1 - r[16*n+8] for n in range(N//32)))
print("2. Theorem 3.9 relations:", t39); ok &= t39
bal = all(sum(r[16*q:16*q+16]) == 8 for q in range(N//16))
print("3. block balance (8 ones per aligned 16-block):", bal); ok &= bal
table = {}
consistent = True
for q in range(N//16):
    key = (r[16*q], q % 4)
    blk = tuple(r[16*q:16*q+16])
    if key in table and table[key] != blk: consistent = False
    table[key] = blk
print("4. block determined by (r(16q), q mod 4):", consistent, "keys:", sorted(table))
ok &= consistent
F = {i: {key: [t for t in range(16) if blk[t] == i] for key, blk in table.items()} for i in (0, 1)}
u = {i: [m for m in range(N) if r[m] == i] for i in (0, 1)}
form = all(u[i][8*q + j] == 16*q + F[i][(r[16*q], q % 4)][j] for i in (0, 1) for q in range(N//16 - 1) for j in range(8))
print("   u_i(8q+j) = 16q + F_i(r(16q), q mod 4, j):", form); ok &= form
for i in (0, 1):
    print(f"   F_{i}:", {k: v for k, v in sorted(F[i].items())})
ranks = {}
for i in (0, 1):
    for E in (6, 7, 8):
        rows = [[u[i][(2**e)*n + b] for n in range(200)] for e in range(E + 1) for b in range(2**e)]
        ranks[(i, E)] = int(np.linalg.matrix_rank(np.array(rows, dtype=float)))
stab = all(ranks[(i, 6)] == ranks[(i, 7)] == ranks[(i, 8)] for i in (0, 1))
print("5. 2-kernel span ranks (E=6,7,8):", ranks, "stable:", stab); ok &= stab
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
