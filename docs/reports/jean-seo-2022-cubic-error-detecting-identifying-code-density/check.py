"""Jean-Seo, "Optimal Error-detection system for Identifying Codes", arXiv:2208.06052v1, Conjecture 1:
for a cubic graph G, DET:IC%(G) >= 3/4.

Search-seat family G_t (t >= 1): base F_t = C_{4t} x K_2 (prism, cubic, 8t vertices, k = 12t edges e_0..e_{k-1}
in sorted order). Vertices A_v (v in F_t), B_j (0 <= j < 6t), P_i, Q_i, R_i (0 <= i < k). Edges: A_u P_i, A_v P_i for
e_i = {u, v}; P_i Q_i; R_i Q_i; R_i Q_{i+1 mod k}; R_i B_{floor(i/2)}; B_{2j} B_{2j+1} (0 <= j < 3t).
Code S_t = all P, Q, R (36t of 50t vertices).

Checks, exact, for t = 1..4: G_t is simple, connected and 3-regular; S_t satisfies Theorem 2.1 (every |N_S[v]| >= 2,
every pair has |N_S[u] - N_S[v]| >= 2 or |N_S[v] - N_S[u]| >= 2); 4|S_t| < 3|V(G_t)|.
Also: the detector-neighbourhood sizes by type (A,B,P,Q,R) = (3,2,2,4,3), and the explicit t = 1 edge list on 0..49.
"""
import itertools
import networkx as nx

def build(t):
    n4 = 4 * t
    base = nx.Graph()
    for a in range(n4):
        for b in (0, 1):
            base.add_edge((a, b), ((a + 1) % n4, b))
        base.add_edge((a, 0), (a, 1))
    E = sorted(tuple(sorted(e)) for e in base.edges())
    k = len(E)
    G = nx.Graph()
    for i, (u, v) in enumerate(E):
        G.add_edge(('A', u), ('P', i)); G.add_edge(('A', v), ('P', i))
        G.add_edge(('P', i), ('Q', i)); G.add_edge(('R', i), ('Q', i))
        G.add_edge(('R', i), ('Q', (i + 1) % k)); G.add_edge(('R', i), ('B', i // 2))
    for j in range(3 * t):
        G.add_edge(('B', 2 * j), ('B', 2 * j + 1))
    S = {v for v in G if v[0] in 'PQR'}
    return G, S, base, k

def detic(G, S):
    N = {v: (set(G[v]) | {v}) & S for v in G}
    dom = all(len(N[v]) >= 2 for v in G)
    pairs = 0; bad = 0
    for u, v in itertools.combinations(sorted(G), 2):
        pairs += 1
        if max(len(N[u] - N[v]), len(N[v] - N[u])) < 2:
            bad += 1
    return dom, pairs, bad, N

ok = True
for t in (1, 2, 3, 4):
    G, S, base, k = build(t)
    n = G.number_of_nodes()
    reg = set(d for _, d in G.degree()) == {3}
    con = nx.is_connected(G)
    simple = nx.number_of_selfloops(G) == 0
    baseok = set(d for _, d in base.degree()) == {3} and base.number_of_nodes() == 8 * t and k == 12 * t
    dom, pairs, bad, N = detic(G, S)
    sizes = {T: sorted(set(len(N[v]) for v in G if v[0] == T)) for T in 'ABPQR'}
    print(f"t={t}: n={n} m={G.number_of_edges()} cubic={reg} connected={con} simple={simple} base_ok={baseok} "
          f"|S|={len(S)} 2-dominated={dom} pairs={pairs} undistinguished={bad} sizes={sizes} 4|S|<3n={4*len(S) < 3*n}")
    ok &= reg and con and simple and baseok and dom and bad == 0 and 4 * len(S) < 3 * n and len(S) * 50 == 36 * n
G, S, base, k = build(1)
order = sorted(G, key=lambda v: ('ABPQR'.index(v[0]), v[1]))
ix = {v: i for i, v in enumerate(order)}
edges = sorted(tuple(sorted((ix[u], ix[v]))) for u, v in G.edges())
print("t=1 vertex order:", [f"{v[0]}{v[1]}" for v in order])
print("t=1 edges (0..49):", edges)
print("t=1 code S:", sorted(ix[v] for v in S))
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
