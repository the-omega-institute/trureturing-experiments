#!/usr/bin/env python3
"""Checks for C. Badica and A. Badica, arXiv:2608.18194v1, Section 8 (second open problem).

d(A, B) = 1 - f(A meet B)/f(A join B) (17). Fidelity: Theorem 1 on the Boolean lattice 2^3 with the modular valuation
|S| + 1, and Theorem 4 (necessity) on N5 and M3. Claim: random strictly positive, monotone, supermodular and
log-submodular valuations on the non-distributive lattices N5 and M3 (rejection sampling; 2^3 is distributive and
already covered by Theorem 2) satisfy the triangle
inequality. Exits 1 if a check fails.
"""
import itertools, math, random, sys
random.seed(20261023)

def lattice(elems, leq):
    meet = {}; join = {}
    for a in elems:
        for b in elems:
            lb = [c for c in elems if leq(c, a) and leq(c, b)]
            ub = [c for c in elems if leq(a, c) and leq(b, c)]
            meet[a, b] = next(c for c in lb if all(leq(d, c) for d in lb))
            join[a, b] = next(c for c in ub if all(leq(c, d) for d in ub))
    return meet, join

def make(elems, rel):
    leq = lambda p, q: p == q or (p, q) in rel
    return elems, leq, *lattice(elems, leq)

N5 = make(['0', 'a', 'b', 'c', '1'], {('0', 'a'), ('0', 'b'), ('0', 'c'), ('0', '1'), ('a', 'c'), ('a', '1'), ('b', '1'), ('c', '1')})
M3 = make(['0', 'x', 'y', 'z', '1'], {('0', 'x'), ('0', 'y'), ('0', 'z'), ('0', '1'), ('x', '1'), ('y', '1'), ('z', '1')})
B3e = [frozenset(s) for r in range(4) for s in itertools.combinations(range(3), r)]
B3 = (B3e, lambda p, q: p <= q, {(a, b): a & b for a in B3e for b in B3e}, {(a, b): a | b for a in B3e for b in B3e})
EPS = 1e-12

def tri_ok(L, f):
    E, leq, meet, join = L
    d = lambda a, b: 1 - f[meet[a, b]] / f[join[a, b]]
    return all(d(A, C) <= d(A, B) + d(B, C) + EPS for A, B, C in itertools.product(E, repeat=3))

def monotone(L, f):
    E, leq, meet, join = L
    return all(f[a] <= f[b] + EPS for a in E for b in E if leq(a, b))

def supermod(L, f):
    E, leq, meet, join = L
    return all(f[a] + f[b] <= f[join[a, b]] + f[meet[a, b]] + EPS for a in E for b in E)

def logsubmod(L, f):
    E, leq, meet, join = L
    return all(math.log(f[join[a, b]]) + math.log(f[meet[a, b]]) <= math.log(f[a]) + math.log(f[b]) + EPS for a in E for b in E)

ok = True
f = {s: len(s) + 1 for s in B3e}
t1 = tri_ok(B3, f)
print(f"fidelity Theorem 1 (2^3, f = |S| + 1): triangle inequality {t1}")
ok &= t1
for name, L in (('N5', N5), ('M3', M3)):
    E = L[0]; tri = 0; bad = 0
    for _ in range(200000):
        f = {e: random.uniform(0.1, 10) for e in E}
        if tri_ok(L, f):
            tri += 1
            g = {e: 1 / f[e] for e in E}
            bad += not (monotone(L, f) and supermod(L, f) and supermod(L, g))
    print(f"fidelity Theorem 4 ({name}): {tri} valuations satisfy the triangle inequality; {bad} fail the necessary conditions")
    ok &= bad == 0
for name, L in (('N5', N5), ('M3', M3)):
    E = L[0]; adm = 0; viol = 0
    for _ in range(400000):
        f = {e: random.uniform(0.1, 10) for e in E}
        if monotone(L, f) and supermod(L, f) and logsubmod(L, f):
            adm += 1
            viol += not tri_ok(L, f)
    print(f"claim ({name}): {adm} admissible valuations, {viol} triangle violations")
    ok &= viol == 0 and adm > 0
print("ALL_CHECKS_OK" if ok else "CHECK_FAIL")
sys.exit(0 if ok else 1)
