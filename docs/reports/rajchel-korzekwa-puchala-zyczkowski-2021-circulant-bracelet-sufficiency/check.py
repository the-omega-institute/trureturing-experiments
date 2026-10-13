#!/usr/bin/env python3
"""Circulant bracelet matrices that are not unistochastic.

Source question: Rajchel-Mieldzioc, Korzekwa, Puchala, Zyczkowski, "Algebraic and geometric
structures inside the Birkhoff polytope", arXiv:2101.11288; J. Math. Phys. 63, 012202 (2022),
Section V: is the bracelet condition sufficient for unistochasticity of circulant matrices, d > 4?

For a support set S of Z_d, B_S is the d x d circulant matrix with B[i][j] = 1/|S| when
(j - i) mod d lies in S and 0 otherwise.  Every entry of B_S is 0 or 1/|S|, so for a pair of
distinct rows (or columns) the bracelet lengths sqrt(B B') are 1/|S| on the common support and
0 elsewhere.  The bracelet inequality 2 max <= sum then holds exactly when the common support
does not have exactly one element.

Support certificate.  Suppose a unitary U has |U_ij|^2 = B_ij, so U has the support of B.
Choose columns p, q whose supports meet exactly in a row set T.  If a row i outside T meets
every row of T in the same two columns {c, c'}, row orthogonality gives
U[r][c'] = lam * U[r][c] for all r in T with one nonzero lam.  If such links chain p to q, then
U[r][q] = mu * U[r][p] on T with mu nonzero, and the orthogonality of columns p and q
(sum over T only) gives conj(mu) * sum_{r in T} |U[r][p]|^2 = 0, a contradiction.

Prints one line per check and ALL_OK on success; exits 1 on the first failure.
"""
import itertools, sys
from fractions import Fraction

def fail(msg):
    print("FAIL", msg); sys.exit(1)

def matrix(S, d):
    w = Fraction(1, len(S))
    return [[w if (j - i) % d in S else Fraction(0) for j in range(d)] for i in range(d)]

def overlap_counts(S, d):
    return [len(S & {(s + h) % d for s in S}) for h in range(1, d)]

def check_bistochastic_circulant(B, d):
    for i in range(d):
        if sum(B[i]) != 1: return False
        if sum(B[r][i] for r in range(d)) != 1: return False
        for j in range(d):
            if B[i][j] < 0 or B[i][j] != B[(i + 1) % d][(j + 1) % d]: return False
    return True

def check_bracelet(B, d):
    # entries are 0 or w, so sqrt(B[a][k] * B[a][l]) is w or 0: compare counts of nonzero lengths
    for k in range(d):
        for l in range(d):
            if k == l: continue
            col = sum(1 for a in range(d) if B[a][k] and B[a][l])
            row = sum(1 for a in range(d) if B[k][a] and B[l][a])
            if col == 1 or row == 1: return False
    return True

def certificate(S, d):
    rows = [{(i + s) % d for s in S} for i in range(d)]
    cols = [{(c - s) % d for s in S} for c in range(d)]
    p = 0
    for q in range(1, d):
        T = cols[p] & cols[q]
        if not T: continue
        over = [c for c in range(d) if T <= cols[c]]
        adj = {c: set() for c in over}
        links = []
        for i in range(d):
            if i in T: continue
            meets = {frozenset(rows[i] & rows[r]) for r in T}
            if len(meets) == 1:
                (m,) = meets
                if len(m) == 2 and all(c in adj for c in m):
                    a, b = sorted(m); adj[a].add(b); adj[b].add(a); links.append((i, a, b))
        seen, stack = {p}, [p]
        while stack:
            u = stack.pop()
            for v in adj[u]:
                if v not in seen: seen.add(v); stack.append(v)
        if q in seen:
            return q, sorted(T), links
    return None

def counterexample(S, d, label):
    B = matrix(S, d)
    if not check_bistochastic_circulant(B, d): fail(f"{label} not circulant bistochastic")
    if not check_bracelet(B, d): fail(f"{label} bracelet")
    if max(B[0]) > Fraction(1, 2): fail(f"{label} entry above 1/2")
    c = certificate(S, d)
    if not c: fail(f"{label} no certificate")
    print(f"{label}: d={d} S={sorted(S)} overlaps={overlap_counts(S, d)} columns=(0,{c[0]}) T={c[1]} links={c[2]}")

# 1. The d = 12 counterexample.
S12 = {0, 1, 3, 4, 6, 9}
counterexample(S12, 12, "D12")
rows = [{(i + s) % 12 for s in S12} for i in range(12)]
cols = [{(c - s) % 12 for s in S12} for c in range(12)]
assert cols[0] & cols[6] == {0, 3, 6, 9}
assert all(rows[11] & rows[r] == {0, 3} for r in (0, 3, 6, 9))
assert all(rows[2] & rows[r] == {3, 6} for r in (0, 3, 6, 9))
print("D12_CHAIN: rows 11 and 2 meet rows 0,3,6,9 in columns {0,3} and {3,6}; columns 0 and 6 meet exactly in rows {0,3,6,9}")

# 2. The d = 15 pattern obtained from the Lundgren-Reid-Severini-Stewart block matrix by the CRT relabelling.
S15 = {1, 4, 6, 7, 10, 12, 13}
crt = {d for d in range(15) if (d % 3 == 0 and d % 5 in (1, 2)) or d % 3 == 1}
if crt != S15: fail("CRT relabelling")
counterexample(S15, 15, "D15")

# 3. Multiples: S scaled by m inside Z_{12m}, Z_{14m}, Z_{15m} (a direct sum of m copies).
for base_d, S in ((12, S12), (14, {0, 1, 2, 3, 6, 7}), (15, S15)):
    for m in (2, 3):
        counterexample({m * s for s in S}, m * base_d, f"D{base_d}x{m}")

# 4. Certificate search over all circulant supports containing 0, d = 3..15.
for d in range(3, 16):
    found = 0
    for k in range(2, d):
        for T in itertools.combinations(range(1, d), k - 1):
            S = {0, *T}
            if 1 in overlap_counts(S, d): continue
            if certificate(S, d): found += 1
    print(f"SEARCH d={d} supports_with_certificate={found}")

print("ALL_OK")
