#!/usr/bin/env python3
"""Bajnok-Berson-Just, arXiv:2211.13675v1 (Involve 15 (2022)), Question 6: is there a noncyclic finite abelian group with a
perfect restricted s-basis of size s + 1 (every element is the sum of exactly one subset T of A with |T| <= s)?
Then |G| = 2^(s+1) - 1. Exhaustive search (pure Python, DFS with distinct-subset-sum pruning) over every abelian group of
order 2^m - 1 for m = 2..6: the cyclic groups Z3, Z7, Z15, Z31, Z63 and the only noncyclic one in range, Z3 x Z21.
Checks: no basis exists in Z3 x Z21; in each cyclic group every basis element is a generator (the parity argument's
conclusion); and counts of bases. usage: check.py"""
import sys, math

def group(mods):
    n = 1
    for q in mods: n *= q
    elems = []
    for i in range(n):
        v, x = [], i
        for q in mods: v.append(x % q); x //= q
        elems.append(tuple(v))
    return elems

def add(a, b, mods): return tuple((x + y) % q for x, y, q in zip(a, b, mods))

def order(a, mods):
    o = 1
    for x, q in zip(a, mods): o = o * (q // math.gcd(x, q)) // math.gcd(o, q // math.gcd(x, q))
    return o

def search(mods, m):
    """All m-subsets A of G (as sorted index tuples) whose proper-subset sums are pairwise distinct (hence a bijection onto G,
    since there are 2^m - 1 = |G| proper subsets)."""
    G = group(mods); idx = {g: i for i, g in enumerate(G)}; zero = tuple(0 for _ in mods)
    found = []
    # sums[k] = list of (sum, size) for all subsets of the chosen prefix; keep only sizes <= m - 1
    def dfs(start, chosen, sums):
        if len(chosen) == m:
            found.append(tuple(chosen)); return
        for i in range(start, len(G)):
            g = G[i]
            if g == zero: continue
            new = [(add(s, g, mods), k + 1) for s, k in sums if k + 1 <= m - 1]
            seen = {s for s, _ in sums}
            ok = True
            for s, _ in new:
                if s in seen: ok = False; break
                seen.add(s)
            if not ok: continue
            dfs(i + 1, chosen + [i], sums + new)
    dfs(0, [], [(zero, 0)])
    return G, found

ok = True
for mods, m in [((3,), 2), ((7,), 3), ((15,), 4), ((31,), 5), ((63,), 6), ((3, 21), 6)]:
    G, found = search(mods, m)
    n = len(G)
    cyclic = len(mods) == 1
    gens_ok = all(order(G[i], mods) == n for A in found for i in A)
    good = (len(found) > 0 and gens_ok) if cyclic else (len(found) == 0)
    ok &= good
    print(f"G = {' x '.join('Z%d' % q for q in mods)} (|G| = {n}, s = {m - 1}): perfect restricted bases of size s+1: {len(found)}; "
          f"every element a generator: {gens_ok}; {'PASS' if good else 'FAIL'}", flush=True)
print("RESULT PASS" if ok else "RESULT FAIL"); print("EXIT=0" if ok else "EXIT=1"); sys.exit(0 if ok else 1)
