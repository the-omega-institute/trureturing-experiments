"""Check for Question 3.17 (v2; 3.16 in v1) of arXiv:2502.10521: the graphs K2 x K5, C4 x C3 (Cartesian products) and
the tesseract Q4.

For each graph: the vertex count, regularity, the rank of the adjacency matrix (the source notes the matrices are
singular), the absence of ordinary twins (distinct vertices with equal open neighbourhoods), and Lemma 1 of the
proof: with h_v the indicator of N(v), the product p_x = prod_{v: x in N(v)} h_v * prod_{v: x notin N(v)} (1 - h_v)
equals the point indicator delta_x for every vertex x.
"""
import itertools
import sympy as sp


def cart(V1, adj1, V2, adj2):
    V = [(a, b) for a in V1 for b in V2]
    adj = lambda u, w: (u[0] == w[0] and adj2(u[1], w[1])) or (u[1] == w[1] and adj1(u[0], w[0]))
    return V, adj


Kn = lambda n: (list(range(n)), lambda i, j: i != j)
Cn = lambda n: (list(range(n)), lambda i, j: (i - j) % n in (1, n - 1))
graphs = {
    'K2 x K5': cart(*Kn(2), *Kn(5)),
    'C4 x C3': cart(*Cn(4), *Cn(3)),
    'Q4': (list(itertools.product((0, 1), repeat=4)), lambda u, w: sum(a != b for a, b in zip(u, w)) == 1),
}
ok = True
for name, (V, adj) in graphs.items():
    n = len(V)
    A = sp.Matrix(n, n, lambda i, j: 1 if adj(V[i], V[j]) else 0)
    assert A == A.T and all(A[i, i] == 0 for i in range(n))
    deg = {sum(A.row(i)) for i in range(n)}
    N = [frozenset(j for j in range(n) if A[i, j]) for i in range(n)]
    twins = [(i, j) for i in range(n) for j in range(i + 1, n) if N[i] == N[j]]
    lemma1 = True
    for x in range(n):
        for y in range(n):
            val = 1
            for v in range(n):
                hv = 1 if y in N[v] else 0
                val *= hv if x in N[v] else (1 - hv)
            lemma1 &= (val == (1 if x == y else 0))
    print(f"{name}: {n} vertices, degrees {sorted(deg)}, rank A = {A.rank()} (singular: {A.rank() < n}), "
          f"ordinary twin pairs = {len(twins)}, Lemma 1 p_x = delta_x for all x: {lemma1}")
    ok &= (len(twins) == 0) and lemma1 and A.rank() < n
# control: C4 = Q2 has twins (opposite vertices), so Lemma 1 must fail there
V, adj = Cn(4)
N = [frozenset(j for j in V if adj(i, j)) for i in V]
print(f"control C4: twin pairs {[(i, j) for i in V for j in V if i < j and N[i] == N[j]]}")
print("ALL_CHECKS", ok)
raise SystemExit(0 if ok else 1)
