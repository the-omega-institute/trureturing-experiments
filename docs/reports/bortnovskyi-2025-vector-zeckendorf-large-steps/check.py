#!/usr/bin/env python3
"""Checks for Bortnovskyi et al., arXiv:2510.15923v2, Conjecture 1.19 (large steps bound j_lsb).

Implements Definition 1.4 (k-bonacci vectors, via X_{-n-k} = X_{-n} - sum_{i=1}^{k-1} X_{-n-i}), Definition 1.15
(small steps), Definition 1.17 (large steps; every minimizing index is recorded, ties are counted) and
j_lsb = k(N - 1) + max n_i. Exits 1 if a fidelity check fails.
"""
import math, sys

def kbonacci(k, N):
    Y = [tuple(0 for _ in range(k - 1))] + [tuple(1 if j == i else 0 for j in range(k - 1)) for i in range(k - 1)]
    while len(Y) < N:
        n = len(Y) - k
        Y.append(tuple(Y[n][c] - sum(Y[n + j][c] for j in range(1, k)) for c in range(k - 1)))
    return Y  # Y[n] = X_{-n}

def n2(v):
    return sum(x * x for x in v)

def small_steps(k, v):
    if any(x < 0 for x in v):
        vm = min(v)
        return [k] * (-vm) + [j + 1 for j in range(k - 1) for _ in range(v[j] - vm)]
    return [j + 1 for j in range(k - 1) for _ in range(v[j])]

def large_steps(k, v, Y):
    ns, ties = [], 0
    while True:
        cur = n2(v)
        cands = [(n2(tuple(a - b for a, b in zip(v, Y[n]))), n) for n in range(1, len(Y))]
        best = min(d for d, _ in cands)
        arg = [n for d, n in cands if d == best]
        if best >= cur:
            return ns, v, ties
        ties += len(arg) > 1
        n = arg[0]
        ns.append(n)
        v = tuple(a - b for a, b in zip(v, Y[n]))

def jlsb(k, ns, res):
    allidx = ns + small_steps(k, res)
    return k * (len(allidx) - 1) + max(allidx, default=0), len(allidx)

def table(k, Y, m):
    while n2(Y[-1]) < 16 * n2(m):
        Y.extend(kbonacci(k, len(Y) + 64)[len(Y):])
    return Y

ok = True
Y3 = kbonacci(3, 200)
assert Y3[3] == (-1, -1) and Y3[4] == (2, 0)
for v, want in (((2, -2), 16), ((3, 0), 7)):
    ns, res, _ = large_steps(3, v, Y3)
    j, N = jlsb(3, ns, res)
    print(f"fidelity k=3 v={v}: accepted={ns} residual={res} small={small_steps(3, res)} jlsb={j} (source {want})")
    ok &= j == want
for k, shape in ((3, lambda m: (0, m)), (4, lambda m: (0, m, 0))):
    Y = kbonacci(k, 60)
    print(f"k={k} v=" + ("(0,m)" if k == 3 else "(0,m,0)"))
    for e in range(4, 41, 4):
        m = 2 ** e
        Y = table(k, Y, shape(m))
        ns, res, ties = large_steps(k, shape(m), Y)
        j, N = jlsb(k, ns, res)
        print(f"  m=2^{e}: accepted={len(ns)} N={N} jlsb={j} jlsb/log(norm)={j / math.log(m):.2f} ties={ties}")
print("FIDELITY_OK" if ok else "FIDELITY_FAIL")
sys.exit(0 if ok else 1)
