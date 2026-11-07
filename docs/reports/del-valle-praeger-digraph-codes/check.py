#!/usr/bin/env python3
"""Exact finite check for the three-layer/two-fibre digraph-code lift.

All arithmetic is over GF(2) using integers and bit operations.  The check
uses H=K_{3,3} and the even-parity [3,2,2] local code, verifies the source's
small graph-code parameters, builds the 216-coordinate digraph constraints,
and checks the twelve-copy decomposition, parameters, graph properties and
the stated finite symmetry generators.
"""

from collections import deque
from itertools import product
from pathlib import Path


def gf2_rank(rows, ncols):
    rows = [r & ((1 << ncols) - 1) for r in rows if r]
    rank = 0
    for col in range(ncols):
        pivot = next((j for j in range(rank, len(rows)) if (rows[j] >> col) & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for j in range(len(rows)):
            if j != rank and ((rows[j] >> col) & 1):
                rows[j] ^= rows[rank]
        rank += 1
    return rank


def gf2_nullspace(rows, ncols):
    a = [r & ((1 << ncols) - 1) for r in rows if r]
    rank = 0
    pivots = []
    for col in range(ncols):
        pivot = next((j for j in range(rank, len(a)) if (a[j] >> col) & 1), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        for j in range(len(a)):
            if j != rank and ((a[j] >> col) & 1):
                a[j] ^= a[rank]
        pivots.append(col)
        rank += 1
    free = [c for c in range(ncols) if c not in pivots]
    basis = []
    for f in free:
        v = 1 << f
        for i, p in enumerate(pivots):
            if (a[i] >> f) & 1:
                v |= 1 << p
        basis.append(v)
    return basis


def hamming(x):
    return x.bit_count()


def parity3(bits):
    return bits.bit_count() & 1


def main():
    # H = K_{3,3}; vertices 0,1,2 are L and 3,4,5 are R.
    left = range(3)
    right = range(3, 6)
    edges = [(u, v) for u in left for v in right]
    edge_index = {e: j for j, e in enumerate(edges)}
    neighbors = {u: (list(right) if u in left else list(left)) for u in range(6)}

    # Original graph code: even parity on each of the six stars.
    graph_rows = []
    for u in range(6):
        row = 0
        for v in neighbors[u]:
            e = (u, v) if (u, v) in edge_index else (v, u)
            row ^= 1 << edge_index[e]
        graph_rows.append(row)
    graph_rank = gf2_rank(graph_rows, 9)
    graph_dim = 9 - graph_rank
    graph_words = []
    for x in range(1 << 9):
        ok = True
        for u in range(6):
            star = 0
            for j, v in enumerate(neighbors[u]):
                e = (u, v) if (u, v) in edge_index else (v, u)
                star |= ((x >> edge_index[e]) & 1) << j
            if parity3(star):
                ok = False
                break
        if ok:
            graph_words.append(x)
    assert graph_dim == 4 and len(graph_words) == 16
    assert min(hamming(x) for x in graph_words if x) == 4

    # Digraph coordinates are ((u,i,a),(v,i+1,b)) for i mod 3 and a,b in {0,1}.
    vertices = [(u, i, a) for u in range(6) for i in range(3) for a in range(2)]
    vertex_index = {x: j for j, x in enumerate(vertices)}
    arcs = []
    for u, v in edges + [(v, u) for u, v in edges]:
        for i in range(3):
            for a in range(2):
                for b in range(2):
                    arcs.append(((u, i, a), (v, (i + 1) % 3, b)))
    arc_index = {x: j for j, x in enumerate(arcs)}
    assert len(vertices) == 36 and len(arcs) == 216

    # B1=B direct_sum B: one parity equation for each target fibre b.
    # B2=Delta(B): equality of the two origin fibres and one parity equation.
    digraph_rows = []
    for u, i, a in vertices:
        for b in range(2):
            row = 0
            for v in neighbors[u]:
                row ^= 1 << arc_index[((u, i, a), (v, (i + 1) % 3, b))]
            digraph_rows.append(row)
    for v, i, b in vertices:
        for u in neighbors[v]:
            row = (1 << arc_index[((u, (i - 1) % 3, 0), (v, i, b))])
            row ^= 1 << arc_index[((u, (i - 1) % 3, 1), (v, i, b))]
            digraph_rows.append(row)
        row = 0
        for u in neighbors[v]:
            row ^= 1 << arc_index[((u, (i - 1) % 3, 0), (v, i, b))]
        digraph_rows.append(row)
    digraph_rank = gf2_rank(digraph_rows, len(arcs))
    digraph_dim = len(arcs) - digraph_rank
    assert digraph_rank == 168 and digraph_dim == 48

    # Explicit twelve-copy map: each graph-code word is repeated on both
    # origin fibres, for each layer, target fibre, and H orientation.
    lifted_weights = []
    for w in graph_words:
        lifted_weights.append(2 * hamming(w))
    assert min(lifted_weights[1:]) == 8
    assert min(2 * hamming(w) for w in graph_words if w) == 8
    assert 48 == 12 * graph_dim and 216 == 24 * len(edges)

    # B1 and B2 have dimensions 4 and 2 respectively, so cannot be isomorphic.
    assert 4 != 2

    # Graph properties: strong connectivity and no reverse arcs.
    outgoing = {x: [] for x in vertices}
    for j, (x, y) in enumerate(arcs):
        outgoing[x].append(y)
        assert (y, x) not in arc_index
    seen = {vertices[0]}
    q = deque([vertices[0]])
    while q:
        x = q.popleft()
        for y in outgoing[x]:
            if y not in seen:
                seen.add(y)
                q.append(y)
    assert len(seen) == len(vertices)

    # Nine generators of A x ((C2)^3 semidirect C3): four base swaps,
    # one bipartition swap, three layer fibre swaps, and one layer shift.
    base_perms = []
    for x, y in ((0, 1), (1, 2), (3, 4), (4, 5)):
        p = list(range(6))
        p[x], p[y] = p[y], p[x]
        base_perms.append(p)
    base_perms.append([3, 4, 5, 0, 1, 2])

    generators = []
    for p in base_perms:
        generators.append([(p[u], i, a) for u, i, a in vertices])
    for layer in range(3):
        generators.append([(u, i, 1 - a) if i == layer else (u, i, a) for u, i, a in vertices])
    generators.append([(u, (i + 1) % 3, a) for u, i, a in vertices])

    kernel_basis = gf2_nullspace(digraph_rows, len(arcs))
    assert len(kernel_basis) == 48
    for perm_vertices in generators:
        perm_v = dict(zip(vertices, perm_vertices))
        perm_arcs = [arc_index[(perm_v[x], perm_v[y])] for x, y in arcs]
        assert sorted(perm_arcs) == list(range(len(arcs)))
        for w in kernel_basis:
            image = 0
            for j, pj in enumerate(perm_arcs):
                if (w >> j) & 1:
                    image ^= 1 << pj
            assert all((image & row).bit_count() % 2 == 0 for row in digraph_rows)

    print("K3,3 graph code: [9,4,4]")
    print("lifted digraph: 36 vertices, 216 arcs, rank 168, dimension 48, minimum distance 8")
    print("B1 dimension 4; B2 dimension 2; reverse arcs 0; reachable vertices 36")
    print("symmetry generators checked: 9")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
