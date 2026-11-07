#!/usr/bin/env python3
"""Exact checks for Fomin, arXiv:2605.05506v3, Question 4.1.

All geometric decisions use integers or fractions.Fraction. No numerical
approximation, external solver, or third-party program is used.
"""

from fractions import Fraction as Q
from itertools import combinations
import json
from pathlib import Path


WITNESS = [
    (9817, 51), (-8208, 6027), (3083, -9321), (3196, 9669),
    (-7912, -5812), (10183, -51), (-7972, 5729), (3098, -9701),
    (2985, 9352), (-8269, -5944),
]

# Literal vertex order in the TeX picture illustrating K10 (k = 2).
# The source states that this figure has type <12 x 4>.
SOURCE_K10 = [
    (4, 3), (5, 8), (6, 1), (5, 7), (3, 2), (5, 6),
    (7, 1), (5, 5), (2, 2), (5, 4), (8, 1), (1, 2),
]


def sub(a, b):
    return a[0] - b[0], a[1] - b[1]


def det(a, b):
    return a[0] * b[1] - a[1] * b[0]


def orient(a, b, c):
    return det(sub(b, a), sub(c, a))


def on_segment(x, a, b):
    return (orient(a, b, x) == 0
            and min(a[0], b[0]) <= x[0] <= max(a[0], b[0])
            and min(a[1], b[1]) <= x[1] <= max(a[1], b[1]))


def crossing(a, b, c, d):
    """Return exact (t,u,X) for a+t(b-a)=c+u(d-c), 0<t,u<1."""
    v, w, z = sub(b, a), sub(d, c), sub(c, a)
    denominator = det(v, w)
    via_orientations = (orient(a, b, c) * orient(a, b, d) < 0
                        and orient(c, d, a) * orient(c, d, b) < 0)
    if denominator == 0:
        assert not via_orientations
        return None
    t = Q(det(z, w), denominator)
    u = Q(det(z, v), denominator)
    via_parameters = 0 < t < 1 and 0 < u < 1
    assert via_parameters == via_orientations
    if not via_parameters:
        return None
    x = (a[0] + t * v[0], a[1] + t * v[1])
    assert x == (c[0] + u * w[0], c[1] + u * w[1])
    return t, u, x


def components(n, edge_list):
    """Components on the polyline's own vertex set, preserving edge incidences."""
    adj = [[] for _ in range(n)]
    for a, b in edge_list:
        assert a != b
        adj[a].append(b)
        adj[b].append(a)
    assert all(len(row) == 2 for row in adj)
    remaining, result = set(range(n)), []
    while remaining:
        todo, seen = [min(remaining)], set()
        while todo:
            a = todo.pop()
            if a in seen:
                continue
            seen.add(a)
            todo.extend(adj[a])
        remaining -= seen
        result.append(sorted(seen))
    return result


def surgery(n, i, j, near_i, near_j):
    """K9/2 literally: remove AB,MN; add AN,BM.

    near_i (resp. near_j) is the endpoint B (resp. N) adjacent to X
    with no intervening crossing. A and M are the other endpoints.
    """
    ends_i, ends_j = (i, (i + 1) % n), (j, (j + 1) % n)
    assert near_i in ends_i and near_j in ends_j
    b, nn = near_i, near_j
    a = next(v for v in ends_i if v != b)
    m = next(v for v in ends_j if v != nn)
    edges = [(h, (h + 1) % n) for h in range(n) if h not in (i, j)]
    edges.extend([(a, nn), (b, m)])
    return {"removed": [list(ends_i), list(ends_j)],
            "added": [[a, nn], [b, m]],
            "components": components(n, edges)}


def analyze(vertices, k, strict_general_position):
    p = [tuple(map(Q, v)) for v in vertices]
    n = len(p)
    assert n >= 3 and len(set(p)) == n
    edges = [(p[i], p[(i + 1) % n]) for i in range(n)]
    assert all(a != b for a, b in edges)
    assert len(components(n, [(i, (i + 1) % n) for i in range(n)])) == 1
    assert edges[-1][1] == edges[0][0]
    vertex_incidences = []
    for i, (a, b) in enumerate(edges):
        for h, x in enumerate(p):
            if h not in (i, (i + 1) % n) and on_segment(x, a, b):
                vertex_incidences.append([i, h])
    assert not vertex_incidences
    rows, points, crossing_data, parallel = [[] for _ in p], {}, [], []
    for i, j in combinations(range(n), 2):
        a, b = edges[i]
        c, d = edges[j]
        if det(sub(b, a), sub(d, c)) == 0:
            parallel.append([i, j])
        hit = crossing(a, b, c, d)
        if hit is None:
            continue
        t, u, x = hit
        assert len({i, (i + 1) % n, j, (j + 1) % n}) == 4
        rows[i].append((t, j))
        rows[j].append((u, i))
        points.setdefault(x, []).append([i, j])
        crossing_data.append({"edges": [i, j], "t": str(t), "u": str(u),
                              "point": list(map(str, x))})
    orders = [[j for _, j in sorted(row)] for row in rows]
    assert [len(row) for row in rows] == [k] * n
    assert len(crossing_data) * 2 == n * k
    if strict_general_position:
        assert not parallel
        assert len(points) == len(crossing_data)
        for a, b, c in combinations(p, 3):
            assert orient(a, b, c) != 0
    candidates, all_crossings = [], []
    for item in crossing_data:
        i, j = item["edges"]
        near_i = ([i] if orders[i][0] == j else [])
        near_i += ([(i + 1) % n] if orders[i][-1] == j else [])
        near_j = ([j] if orders[j][0] == i else [])
        near_j += ([(j + 1) % n] if orders[j][-1] == i else [])
        k1 = bool(near_i and near_j)
        all_crossings.append({"edges": [i, j], "K9_1": k1})
        for b in near_i:
            for nn in near_j:
                s = surgery(n, i, j, b, nn)
                # Independent check of the source's opposite-orientation criterion.
                opposite = ((b == (i + 1) % n) != (nn == (j + 1) % n))
                assert (len(s["components"]) == 1) == opposite
                candidates.append({"edges": [i, j], "near_endpoints": [b, nn],
                                   "positions": ["first" if b == i else "last",
                                                 "first" if nn == j else "last"],
                                   "K9_2": len(s["components"]) == 1, **s})
    return {"n": n, "k": k, "closed_single_cycle": True,
            "vertices": [list(v) for v in vertices],
            "crossings_per_edge": [len(row) for row in rows],
            "proper_crossing_pairs": len(crossing_data),
            "distinct_crossing_points": len(points),
            "nonincident_vertex_on_edge": vertex_incidences,
            "parallel_edge_pairs": parallel,
            "crossing_order": orders,
            "crossing_parameters": [[str(t) for t, _ in sorted(row)] for row in rows],
            "crossings": crossing_data, "K9_1_screen": all_crossings,
            "K9_1_candidates_and_surgeries": candidates,
            "crossings_satisfying_both": sum(c["K9_2"] for c in candidates)}


def main():
    fidelity = analyze(SOURCE_K10, 4, False)
    assert fidelity["proper_crossing_pairs"] == 24
    witness = analyze(WITNESS, 5, True)
    candidates = witness["K9_1_candidates_and_surgeries"]
    assert [c["edges"] for c in candidates] == [[i, i + 5] for i in range(5)]
    assert all(c["positions"] == ["first", "first"] for c in candidates)
    assert all(sorted(map(len, c["components"])) == [5, 5] for c in candidates)
    assert witness["crossings_satisfying_both"] == 0
    # Positive control: K9/2 must also accept a connected opposite-orientation
    # switch. This checks the actual graph operation, without requiring K9/1.
    assert len(surgery(10, 0, 5, 1, 5)["components"]) == 1
    result = {"arxiv": "2605.05506v3", "question": "4.1",
              "arithmetic": "integers and fractions.Fraction; exact",
              "source_fidelity_K10": fidelity,
              "counterexample": witness,
              "family_status": "mechanism proposed; not certified by this finite check"}
    encoded = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    Path(__file__).with_name("result.json").write_text(encoded)
    print("FIDELITY K10: n=12, k=4, proper crossings=24, distinct points=24")
    print("WITNESS: n=10, k=5, proper crossings=25, distinct points=25")
    print("Distinct vertices; no collinear triples, parallel edges, or vertex on another edge")
    for i, order in enumerate(witness["crossing_order"]):
        print(f"EDGE {i}: {order}")
    for c in candidates:
        print(f"K9/1 {c['edges']}: {c['positions']}; "
              f"K9/2={c['K9_2']}; cycles={c['components']}")
    print("K9/2 positive control: connected; K10 example: 7 successful candidates")
    print("Full exact crossing parameters and surgery incidences: result.json")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
