#!/usr/bin/env python3
"""Apte, Parekh, Sud, "Conjectured Bounds for 2-Local Hamiltonians via Token Graphs", arXiv:2506.03441v2,
Appendix (Upper bounds based on maximum cut), Conjecture 7 (conj:A_cut_bound): for all G and 1 <= k <= floor(n/2),
lambda_min(A(F_k(G))) >= -(C(G) + M_k(G))/2, with C the maximum cut weight and M_k the maximum weight of a matching
with at most k edges.

Family: G_s = C_5 blow-up (each cycle vertex -> independent set of size s, each cycle edge -> K_{s,s}), N = 5s,
k = N/2 (s even). Test vector w(S) = prod_{v in S} z^{block(v)}, z = (-4+3i)/5.

Checks:
  1. The Rayleigh-quotient formula R_s = -2527 s^3 / (250 (5s - 1)) against the brute-force quotient on the token
     graph for s = 2 (N = 10, k = 5, 252 vertices) and s = 4 (N = 20, k = 10, 184756 vertices; sparse).
  2. The identity 4 Re z + Re z^4 = -2527/625 (exact, Fractions).
  3. Max cut C(G_s) = 4 s^2: brute force over all 2^10 cuts for s = 2, and the multilinear vertex argument for every s
     (all 32 block patterns of C_5 give at most 4 cut block-edges).
  4. M_k(G_s) <= k (unit weights) and = 5s/2 (explicit perfect matching).
  5. The violation: R_s < -(4 s^2 + 5 s/2)/2 exactly when 108 s^2 - 4250 s + 1250 > 0; true for every even s >= 40 and
     false for s <= 38 (exact rational arithmetic up to s = 200); s = 40 gives R = -646912/199 < -3250.
  6. lambda_min(A(F_5(G_2))) by dense diagonalization (the conjecture holds there, consistent with the source's
     verification up to order ten)."""
import sys, itertools, math
from fractions import Fraction as F
import numpy as np

bad = 0
def report(label, ok):
    global bad
    print(label, ok)
    bad += 0 if ok else 1

def blowup(s):
    n = 5 * s; block = [v // s for v in range(n)]
    edges = [(u, v) for u in range(n) for v in range(u + 1, n) if (block[u] - block[v]) % 5 in (1, 4)]
    return n, block, edges

z = complex(-4, 3) / 5
def R_formula(s): return F(-2527 * s**3, 250 * (5 * s - 1))

# 2. exact identity: z = cos t + i sin t with cos t = -4/5
c = F(-4, 5); cos4 = 8 * c**4 - 8 * c**2 + 1
report('2. 4 Re z + Re z^4 = -2527/625:', 4 * c + cos4 == F(-2527, 625))

# 1. brute-force Rayleigh quotient
def rayleigh(s):
    n, block, edges = blowup(s); k = n // 2
    subsets = list(itertools.combinations(range(n), k)); idx = {S: i for i, S in enumerate(subsets)}
    w = np.array([np.prod([z ** block[v] for v in S]) for S in subsets])
    adj = {v: [] for v in range(n)}
    for u, v in edges: adj[u].append(v); adj[v].append(u)
    num = 0j
    for S in subsets:
        Sset = set(S); iS = idx[S]
        for u in S:
            for v in adj[u]:
                if v in Sset: continue
                T = tuple(sorted((Sset - {u}) | {v}))
                num += np.conj(w[iS]) * w[idx[T]]
    return num / np.vdot(w, w), len(subsets)
ok = True
for s in (2, 4):
    q, m = rayleigh(s)
    ok &= abs(q.imag) < 1e-7 * abs(q.real) and abs(q.real - float(R_formula(s))) < 1e-7 * abs(q.real)   # float products of z powers
    print(f'   s={s}: token-graph order {m}, Rayleigh quotient {q.real:.9f}, formula {float(R_formula(s)):.9f}')
report('1. brute-force Rayleigh quotient equals -2527 s^3/(250(5s-1)) for s = 2, 4:', ok)

# 3. max cut
n, block, edges = blowup(2)
best = max(sum(1 for u, v in edges if ((mask >> u) & 1) != ((mask >> v) & 1)) for mask in range(1 << n))
pat = max(sum(1 for i in range(5) if p[i] != p[(i + 1) % 5]) for p in itertools.product((0, 1), repeat=5))
report('3. C(G_2) = 16 = 4 s^2 by brute force, and C_5 block patterns cut at most 4 block-edges:', best == 16 and pat == 4)

# 4. matching
ok = True
for s in (2, 4, 6, 40):
    n, block, edges = blowup(s); E = set(edges)
    # pair block 0 with 1, 2 with 3, block 4 with itself impossible (independent) -> use 4-0 and 1-2 alternation
    M = []
    half = s // 2
    for t in range(s):
        M.append((0 * s + t, 1 * s + t)) if t < half else M.append((1 * s + t, 2 * s + t))
    for t in range(s):
        M.append((2 * s + t, 3 * s + t)) if t < half else M.append((3 * s + t, 4 * s + t))
    for t in range(half):
        M.append((0 * s + half + t, 4 * s + t))
    used = [v for e in M for v in e]
    ok &= len(used) == len(set(used)) == n and all((min(a, b), max(a, b)) in E for a, b in M) and len(M) == 5 * s // 2
report('4. explicit perfect matching of G_s with 5s/2 = k edges (s = 2, 4, 6, 40), so M_k = k:', ok)

# 5. threshold
ok = True
for s in range(2, 201, 2):
    viol = R_formula(s) < -(F(4 * s * s) + F(5 * s, 2)) / 2
    ok &= viol == (108 * s * s - 4250 * s + 1250 > 0) and viol == (s >= 40)
ok &= R_formula(40) == F(-646912, 199) and F(-646912, 199) < -3250
report('5. violation iff 108 s^2 - 4250 s + 1250 > 0 iff s >= 40 (even s <= 200); R_40 = -646912/199 < -3250:', ok)

# 6. small-case eigenvalue
n, block, edges = blowup(2); k = 5
subsets = list(itertools.combinations(range(n), k)); idx = {S: i for i, S in enumerate(subsets)}
A = np.zeros((len(subsets), len(subsets)))
adjs = {v: set() for v in range(n)}
for u, v in edges: adjs[u].add(v); adjs[v].add(u)
for S in subsets:
    Sset = set(S)
    for u in S:
        for v in adjs[u] - Sset:
            A[idx[S], idx[tuple(sorted((Sset - {u}) | {v}))]] = 1
lam = np.linalg.eigvalsh(A).min()
print(f'   s=2: lambda_min = {lam:.6f}, bound -(16 + 5)/2 = -10.5')
report('6. at s = 2 the conjectured bound holds (lambda_min >= -10.5):', lam >= -10.5 - 1e-9)
print('bad=', bad)
sys.exit(1 if bad else 0)
