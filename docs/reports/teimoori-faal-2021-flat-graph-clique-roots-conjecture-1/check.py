"""Teimoori Faal, arXiv:2112.09721v1, Conjecture 1: a ring of K4 blocks refutes it.

For r >= 4 the graph R_r has shared vertices s_0..s_{r-1} (indices mod r) and private vertices a_i, b_i, with all six
edges on each block {s_i, s_{i+1}, a_i, b_i}. Checks for r = 4..8, by direct enumeration: R_r is 2-connected (connected
after deleting any vertex, at least 3 vertices), K5-free, flat (every edge in at most two triangles), has no isolated
triangle (every triangle has a common neighbour), and its clique polynomial C(x) = 1 + 3r x + 6r x^2 + 4r x^3 + r x^4
has exactly two real roots, both simple, so two nonreal roots; C(-1) = 1, so -1 is not a root.
"""
import itertools
import sys
import numpy as np


def ring(r):
    adj = {}
    for i in range(r):
        block = [("s", i), ("s", (i + 1) % r), ("a", i), ("b", i)]
        for u, v in itertools.combinations(block, 2):
            adj.setdefault(u, set()).add(v)
            adj.setdefault(v, set()).add(u)
    return adj


def connected(adj, removed=None):
    nodes = [v for v in adj if v != removed]
    seen = {nodes[0]}
    stack = [nodes[0]]
    while stack:
        u = stack.pop()
        for w in adj[u]:
            if w != removed and w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == len(nodes)


def cliques(adj, k):
    return [c for c in itertools.combinations(sorted(adj), k) if all(v in adj[u] for u, v in itertools.combinations(c, 2))]


ok = True
for r in range(4, 9):
    adj = ring(r)
    n = len(adj)
    biconn = n >= 3 and connected(adj) and all(connected(adj, v) for v in adj)
    k5free = not cliques(adj, 5)
    tri = cliques(adj, 3)
    flat = all(len(adj[u] & adj[v]) <= 2 for u in adj for v in adj[u])
    noiso = all(adj[t[0]] & adj[t[1]] & adj[t[2]] for t in tri)
    coeffs = [1] + [len(cliques(adj, k)) for k in range(1, 5)]
    expected = [1, 3 * r, 6 * r, 4 * r, r]
    roots = np.roots(coeffs[::-1])
    nonreal = sum(1 for z in roots if abs(z.imag) > 1e-9)
    at_minus_one = sum(c * (-1) ** j for j, c in enumerate(coeffs))
    good = biconn and k5free and flat and noiso and coeffs == expected and nonreal == 2 and at_minus_one == 1
    ok = ok and good
    print(f"r={r}: n={n}, 2-connected {biconn}, K5-free {k5free}, flat {flat}, no isolated triangle {noiso}, "
          f"C = {coeffs}, nonreal roots {nonreal}, C(-1) = {at_minus_one}")
print("RESULT PASS" if ok else "RESULT FAIL")
sys.exit(0 if ok else 1)
