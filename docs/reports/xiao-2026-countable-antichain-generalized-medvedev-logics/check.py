#!/usr/bin/env python3
"""Exact finite checks for Xiao's countable-antichain question.

All searches exhaust their stated finite domain. No numerical approximation,
external solver, or third-party program source is used.
"""
from itertools import combinations, product


def h_edge(n, a, b):
    return a != b and (a - b) % n not in (1, n - 1)


def incidence(n):
    # Root, unordered graph edges, vertices; order by reverse inclusion of
    # their maximal-vertex supports. n >= 5 makes all these supports distinct.
    assert n >= 5
    supports = [frozenset(range(n))]
    supports += [frozenset((a, b)) for a, b in combinations(range(n), 2)
                 if h_edge(n, a, b)]
    supports += [frozenset((a,)) for a in range(n)]
    return supports, [[a >= b for b in supports] for a in supports]


def rooted_pmaps(src_le, dst_le):
    """All root-preserving p-maps between finite posets indexed bottom first.

For every x, the p-map condition is exactly f(up(x)) = up(f(x)).
Assign strict successors first; test that equality at each assigned point.
The root is sent to the root, so the equality there enforces surjectivity.
"""
    ns, nt = len(src_le), len(dst_le)
    src_up = [{j for j in range(ns) if src_le[i][j]} for i in range(ns)]
    dst_up = [{j for j in range(nt) if dst_le[i][j]} for i in range(nt)]
    assert all(i in src_up[i] for i in range(ns))
    assert all(j > i for i in range(ns) for j in src_up[i] - {i})
    f = [-1] * ns
    solutions, nodes = 0, 0

    def visit(x):
        nonlocal solutions, nodes
        nodes += 1
        if x < 0:
            assert set(f) == set(range(nt))
            assert all({f[z] for z in src_up[i]} == dst_up[f[i]]
                       for i in range(ns))
            solutions += 1
            return
        above = {f[z] for z in src_up[x] - {x}}
        for y in ([0] if x == 0 else range(nt)):
            if dst_up[y] == above | {y}:
                f[x] = y
                visit(x - 1)
        f[x] = -1

    visit(ns - 1)
    return solutions, nodes


def color_obstruction(n, m, group_sizes):
    """Exhaustive necessary-condition CSP for maximal-point colorings.

Each group is one active product coordinate. Every pair from different
coordinates must have equal or H_m-adjacent colors. In a root-coordinate,
H_n-adjacent vertices have the same necessary constraint. All m target
colors must occur. Return existence and the number of search nodes.

Domains are integer bitsets. Arc consistency and the requirement that every
color is used are logically necessary pruning rules. Rotating C_m lets us
fix the first variable to color zero without changing existence.
"""
    groups = [(g, j) for g, size in enumerate(group_sizes) for j in range(size)]
    nv = len(groups)
    full = (1 << m) - 1
    neighbors = [set() for _ in groups]
    for a, b in combinations(range(nv), 2):
        g, v = groups[a]
        h, w = groups[b]
        if g != h or (group_sizes[g] == n and h_edge(n, v, w)):
            neighbors[a].add(b)
            neighbors[b].add(a)
    allowed = [sum(1 << b for b in range(m)
                   if a == b or h_edge(m, a, b)) for a in range(m)]
    nodes = 0

    def bits(mask):
        while mask:
            bit = mask & -mask
            yield bit.bit_length() - 1
            mask -= bit

    def visit(dom):
        nonlocal nodes
        nodes += 1
        changed = True
        while changed:
            changed = False
            if any(d == 0 for d in dom):
                return False
            union = 0
            for d in dom:
                union |= d
            if union != full:
                return False
            for color in range(m):
                spots = [i for i, d in enumerate(dom) if d & (1 << color)]
                if len(spots) == 1 and dom[spots[0]] != (1 << color):
                    dom[spots[0]] = 1 << color
                    changed = True
            for i in range(nv):
                for j in neighbors[i]:
                    reduced = sum(1 << c for c in bits(dom[i])
                                  if allowed[c] & dom[j])
                    if reduced != dom[i]:
                        dom[i] = reduced
                        changed = True
        free = [i for i, d in enumerate(dom) if d & (d - 1)]
        if not free:
            return True
        i = min(free, key=lambda j: (dom[j].bit_count(), -len(neighbors[j])))
        for c in bits(dom[i]):
            nxt = dom.copy()
            nxt[i] = 1 << c
            if visit(nxt):
                return True
        return False

    dom = [full] * nv
    dom[0] = 1
    return visit(dom), nodes


def main():
    for k in range(1, 7):
        nonempty = range(1, 1 << k)
        for x, y in product(nonempty, repeat=2):
            fx = tuple(1 if x & (1 << i) else 2 for i in range(k))
            fy = tuple(1 if y & (1 << i) else 2 for i in range(k))
            assert fx != (2,) * k
            assert ((x & y) == y) == all(a <= b for a, b in zip(fx, fy))
        print(f'FIDELITY Proposition 3.1: P_0({k}) ~= C(2)^({k}-), points={(1 << k) - 1}')
    for n in (5, 6):
        a, le = incidence(n)
        assert all(le[i][i] for i in range(len(a)))
        assert all(not (le[i][j] and le[j][i]) or i == j
                   for i, j in product(range(len(a)), repeat=2))
        assert all(not (le[i][j] and le[j][k]) or le[i][k]
                   for i, j, k in product(range(len(a)), repeat=3))
        assert len(a) + 1 == 2 + n + n * (n - 3) // 2
        print(f'FRAME n={n}: H edges={n * (n - 3) // 2}, A points={len(a)}, F points={len(a) + 1}')
    for n, m in product((5, 6), repeat=2):
        src, src_le = incidence(n)
        dst, dst_le = incidence(m)
        count, nodes = rooted_pmaps(src_le, dst_le)
        assert count == (2 * n if n == m else 0)
        # Every non-root principal upset has at most three points. Thus it
        # cannot surject onto A_m. This covers all generated subframes too:
        # restrict any surjective p-map to a preimage of the target root.
        small_cones = [sum(row) for row in src_le[1:]]
        assert max(small_cones) <= 3 < len(dst)
        print(f'EXHAUSTIVE k=1: generated subframe of F_{n}^(-) ->> A_{m}: '
              f'rooted maps={count}, DFS nodes={nodes}')
        for sizes in ((1, n), (2, n), (n, n)):
            exists, states = color_obstruction(n, m, sizes)
            assert not exists
            print(f'EXHAUSTIVE k=2 necessary coloring: n={n}, m={m}, '
                  f'active coatom groups={sizes}, feasible=0, DFS nodes={states}')
    # The other two-active-coordinate cone types have at most four maxima;
    # no such frame can p-map onto A_m with m >= 5 maximal points.
    print('k=2 other active group types (1,1),(1,2),(2,2): <=4 colors, impossible for m>=5')
    print('Finite results test the proposed obstruction; they do not prove the unbounded-power theorem.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
