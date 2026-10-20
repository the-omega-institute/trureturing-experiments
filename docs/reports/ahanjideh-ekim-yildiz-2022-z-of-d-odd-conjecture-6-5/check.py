"""Ahanjideh-Ekim-Yildiz, "Maximum size of a triangle-free graph with bounded maximum degree and matching number",
arXiv:2207.02271v1, Conjecture 6.5: for odd d >= 21, Z(d) = floor(5(d+1)/4), where Z(d) (Definition 4.2) is the least
matching number of an almost d-regular (one vertex of degree d-1, all others d) triangle-free factor-critical graph.

Checks:
 1. The graph B_d of Figure 2 (reconstructed: C5 order A1-A3-A5-A4-A2, bag sizes as in the figure, removed perfect
    matchings A11-A21, A11-A31, A12-A22, A13-A32, A23-A41, A33-A51, A42-A52, and for d = 4k+3 also A41-A51 and a second
    A42-A52) for d = 21, 23, 25, 27, 29, 31: almost d-regular, triangle-free, factor-critical (every G - v has a perfect
    matching), matching number floor(5(d+1)/4).  (Upper bound Z(d) <= Q.)
 2. The finite case analyses of the search seat's lower-bound proof, exhaustively:
    a. d = 23, N = 59 extension: for every a in Z_{>=0}^5 with sum 6 and every w in Z_{>=0}^5 with sum 22,
       if w_{i-1} + w_{i+1} <= a_i whenever w_i > 0 and w_i <= s_i = 14 - a_{i+2} - a_{i+3}, then supp(w) is an
       independent set of C5 (so v can be placed in a class).  Equivalently no such w has an edge in its support.
    b. U = 1 and U = 3 obstructions: for C in 11..60 (U = 1) and 13..60 (U = 3), every u >= 0 with sum U, every
       placement of epsilon (one class gets 1), D_i = s_i u_i + eps_i with s_i = C - u_{i+2} - u_{i+3}, there is an
       independent set I of C5 with sum_{I} D_i > sum_{not I} D_i (so no missing-edge assignment exists).
    c. Arithmetic: for odd d in 21..401 with N <= 2Q - 1 odd: 8(d-1) > 3N except (d, N) = (23, 59); 3(d-1) - N >= 7;
       2N - 5d >= 0 forces (U, C) = (1, 2k+1) for d = 4k+1 and (3, 2k+3) for d = 4k+3.
"""
import itertools, networkx as nx

def build(d):
    k, r = divmod(d, 4)
    assert r in (1, 3)
    if r == 1:
        sizes = {'A11': 1, 'A12': k, 'A13': k, 'A21': 1, 'A22': k, 'A23': k, 'A31': 1, 'A32': k, 'A33': k,
                 'A41': k, 'A42': k + 1, 'A51': k, 'A52': k + 1}
        removed = [('A11', 'A21'), ('A11', 'A31'), ('A12', 'A22'), ('A13', 'A32'), ('A23', 'A41'), ('A33', 'A51'), ('A42', 'A52')]
    else:
        sizes = {'A11': 1, 'A12': k, 'A13': k, 'A21': 1, 'A22': k, 'A23': k + 1, 'A31': 1, 'A32': k, 'A33': k + 1,
                 'A41': k + 1, 'A42': k + 2, 'A51': k + 1, 'A52': k + 2}
        removed = [('A11', 'A21'), ('A11', 'A31'), ('A12', 'A22'), ('A13', 'A32'), ('A23', 'A41'), ('A33', 'A51'),
                   ('A41', 'A51'), ('A42', 'A52'), ('A42', 'A52', 'shift')]
    cls = {b: int(b[1]) for b in sizes}
    adjc = {(1, 3), (3, 5), (5, 4), (4, 2), (2, 1)}
    adjc |= {(b, a) for a, b in adjc}
    verts = [(b, j) for b in sizes for j in range(sizes[b])]
    G = nx.Graph(); G.add_nodes_from(verts)
    for u, v in itertools.combinations(verts, 2):
        if (cls[u[0]], cls[v[0]]) in adjc:
            G.add_edge(u, v)
    for rm in removed:
        X, Y = rm[0], rm[1]; sh = 1 if len(rm) == 3 else 0
        assert sizes[X] == sizes[Y]
        for j in range(sizes[X]):
            e = ((X, j), (Y, (j + sh) % sizes[Y]))
            assert G.has_edge(*e), (d, rm, j)
            G.remove_edge(*e)
    return G

def factor_critical(G):
    for v in G:
        H = G.copy(); H.remove_node(v)
        if len(nx.max_weight_matching(H, maxcardinality=True)) * 2 != H.number_of_nodes():
            return False
    return True

ok = True
for d in (21, 23, 25, 27, 29, 31):
    G = build(d); Q = 5 * (d + 1) // 4
    degs = sorted(dict(G.degree()).values())
    almost = degs[0] == d - 1 and degs[1] == d and degs[-1] == d
    tri = sum(nx.triangles(G).values()) == 0
    fc = factor_critical(G)
    nu = len(nx.max_weight_matching(G, maxcardinality=True))
    print(f"d={d}: n={G.number_of_nodes()} almost-{d}-regular={almost} triangle-free={tri} factor-critical={fc} nu={nu} Q={Q}")
    ok &= almost and tri and fc and nu == Q and G.number_of_nodes() == 2 * Q + 1

def comps(total, parts):
    if parts == 1:
        yield (total,); return
    for x in range(total + 1):
        for rest in comps(total - x, parts - 1):
            yield (x,) + rest

edges5 = [(i, (i + 1) % 5) for i in range(5)]
# 2a
bad = 0; feasible = 0
for a in comps(6, 5):
    s = [14 - a[(i + 2) % 5] - a[(i + 3) % 5] for i in range(5)]
    for w in comps(22, 5):
        if any(w[i] > s[i] for i in range(5)):
            continue
        if any(w[i] > 0 and w[(i - 1) % 5] + w[(i + 1) % 5] > a[i] for i in range(5)):
            continue
        feasible += 1
        supp = {i for i in range(5) if w[i] > 0}
        if any(x in supp and y in supp for x, y in edges5):
            bad += 1
print(f"2a. d=23,N=59: feasible (a,w) = {feasible}, with an edge in supp(w) = {bad}")
ok &= bad == 0 and feasible > 0
# 2b
indep = [I for r in range(1, 3) for I in itertools.combinations(range(5), r)
         if not any(((x, y) in edges5 or (y, x) in edges5) for x, y in itertools.combinations(I, 2))]
surv = 0
for U, Cmin in ((1, 11), (3, 13)):
    for C in range(Cmin, 61):
        for u in comps(U, 5):
            s = [C - u[(i + 2) % 5] - u[(i + 3) % 5] for i in range(5)]
            for e in range(5):
                D = [s[i] * u[i] + (1 if i == e else 0) for i in range(5)]
                if not any(sum(D[i] for i in I) > sum(D[i] for i in range(5) if i not in I) for I in indep):
                    surv += 1
print(f"2b. U in (1,3): configurations with no violated independent-set inequality = {surv}")
ok &= surv == 0
# 2c
arith = True
for d in range(21, 402, 2):
    k, r = divmod(d, 4); Q = 5 * (d + 1) // 4
    for N in range(1, 2 * Q, 2):
        if 2 * N < 5 * (d - 1):
            continue  # a_i >= 0 forces 2N >= 5(d-1)
        hag = 8 * (d - 1) > 3 * N
        if not hag and (d, N) != (23, 59):
            arith = False
        if 3 * (d - 1) - N < 7:
            arith = False
        U = 2 * N - 5 * d
        if U >= 0:
            if r == 1 and (U, N - 2 * d) != (1, 2 * k + 1): arith = False
            if r == 3 and (U, N - 2 * d) != (3, 2 * k + 3): arith = False
print("2c. arithmetic bounds for odd d in 21..401:", arith)
ok &= arith
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
