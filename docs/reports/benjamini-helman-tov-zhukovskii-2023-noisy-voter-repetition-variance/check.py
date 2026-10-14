#!/usr/bin/env python3
"""Linear variance of the repetition count of the discrete-time noisy voter model.

Source: I. Benjamini, H. Helman Tov, M. Zhukovskii, "Global information from local observations of the noisy voter model
on a graph", arXiv:2207.01224v2, Ann. Probab. 51 (2023), Sec. 2.2: with zeta_i = I(X_i(u) = X_{i+d}(u)),
"May Var(zeta_1 + ... + zeta_t) be sublinear? We do not have an answer."

Claim checked (proved in prose by the search seat): Var(W_t) >= (eps^(D+2)/2) floor(t/(d+1)), D = deg(u).
Exact rational computation of Var(W_t) by forward recursion over the joint law of the last d+1 configurations and the
running sum, for small graphs: K2 (u = 0), path P3 (u = end and u = middle), triangle K3, star K_{1,3} (u = leaf, center),
d in {1, 2}, eps in {1/2, 1/3}, t up to 10 (d = 1) or 7 (d = 2). Prints the ratio Var(W_t)/t and checks the bound.
Prints ALL_OK.
"""
import sys, itertools
from fractions import Fraction as Fr
from math import floor

def fail(m):
    print("FAIL", m); sys.exit(1)

def step_law(adj, eps, conf):
    """distribution of the next configuration (tuple of bits) given the current one"""
    n = len(adj)
    per = []
    for v in range(n):
        p1 = eps * Fr(1, 2)                                   # fresh opinion 1
        nb = adj[v]
        p1 += (1 - eps) * Fr(sum(conf[w] for w in nb), len(nb))
        per.append(p1)
    out = {}
    for bits in itertools.product([0, 1], repeat=n):
        p = Fr(1)
        for v in range(n):
            p *= per[v] if bits[v] else 1 - per[v]
        if p: out[bits] = out.get(bits, 0) + p
    return out

def variance(adj, u, eps, d, t):
    n = len(adj)
    # state: tuple of the last d+1 configurations (oldest first); value: dict W -> prob
    init = {}
    for c0 in itertools.product([0, 1], repeat=n):
        init[(c0,)] = {0: Fr(1, 2 ** n)}
    states = init
    # grow to windows of length d+1 without counting
    for _ in range(d):
        new = {}
        for win, wd in states.items():
            for nxt, p in step_law(adj, eps, win[-1]).items():
                key = win + (nxt,)
                dd = new.setdefault(key, {})
                for w, q in wd.items(): dd[w] = dd.get(w, 0) + q * p
        states = new
    # now window = (X_0..X_d); count zeta_i for i = 1..t needs windows (X_i..X_{i+d})
    for i in range(1, t + 1):
        new = {}
        for win, wd in states.items():
            for nxt, p in step_law(adj, eps, win[-1]).items():
                key = win[1:] + (nxt,)
                z = 1 if key[0][u] == key[-1][u] else 0
                dd = new.setdefault(key, {})
                for w, q in wd.items(): dd[w + z] = dd.get(w + z, 0) + q * p
        states = new
    EW = sum(w * q for wd in states.values() for w, q in wd.items())
    EW2 = sum(w * w * q for wd in states.values() for w, q in wd.items())
    tot = sum(q for wd in states.values() for q in wd.values())
    if tot != 1: fail("normalization")
    return EW2 - EW * EW

graphs = {
    "K2": ([[1], [0]], [0]),
    "P3": ([[1], [0, 2], [1]], [0, 1]),
    "K3": ([[1, 2], [0, 2], [0, 1]], [0]),
    "K13": ([[1, 2, 3], [0], [0], [0]], [0, 1]),
}
for name, (adj, us) in graphs.items():
    for u in us:
        Dg = len(adj[u])
        for eps in (Fr(1, 2), Fr(1, 3)):
            for d, T in ((1, 10), (2, 7)):
                if len(adj) == 4 and d == 2: T = 6
                ratios = []
                for t in range(1, T + 1):
                    V = variance(adj, u, eps, d, t)
                    lb = eps ** (Dg + 2) / 2 * floor(t / (d + 1))
                    if V < lb: fail(f"{name} u={u} eps={eps} d={d} t={t}: Var={V} < {lb}")
                    ratios.append(float(V) / t)
                print(f"{name} u={u} eps={eps} d={d}: Var(W_t)/t for t=1..{T}: " + ", ".join(f"{r:.4f}" for r in ratios))
print("ALL_OK")
