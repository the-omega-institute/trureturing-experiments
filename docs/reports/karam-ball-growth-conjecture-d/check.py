#!/usr/bin/env python3
"""Exact checks for Karam, arXiv:2610.20088v1, Conjecture D.

No numerical eigensolver or floating-point arithmetic is used. Lifted edges are
enumerated by their nonbacktracking edge histories, so distinct lifts are never
identified merely because their projections coincide in the finite graph.
"""
from fractions import Fraction as F

A = (3, -8, 10, -8, 3)
T = (6, -5, 2, 2, -5)
U = (-4, 3, -1, -1, 3)


def model(c, radius):
    assert c > 0 and radius >= 0
    k = radius // c
    return 3 * c * (2**k - 1) + 3 * (radius - k * c) * 2**k


class Graph:
    def __init__(self, vertices, edges):
        self.vertices = tuple(vertices)
        self.edges = tuple(edges)
        self.adj = {v: [] for v in self.vertices}
        for i, (v, w, length) in enumerate(self.edges):
            assert length > 0
            self.adj[v].append((i, w, length))
            self.adj[w].append((i, v, length))

    def ball(self, root, radius):
        """Length of the lifted closed ball: sum full and clipped lifted edges."""
        assert radius >= 0
        total = F(0)
        stack = [(root, None, radius)]
        while stack:
            v, incoming, remaining = stack.pop()
            for edge, w, length in self.adj[v]:
                if edge == incoming:
                    continue
                total += min(length, remaining)
                if length < remaining:
                    stack.append((w, edge, remaining - length))
        return total

    def connected(self):
        reached = {self.vertices[0]}
        stack = list(reached)
        while stack:
            for _, w, _ in self.adj[stack.pop()]:
                if w not in reached:
                    reached.add(w)
                    stack.append(w)
        return len(reached) == len(self.vertices)


def family(d, eps, scale):
    assert d >= 1 and 0 <= eps < F(1, 60) and scale > 0
    vertices = [(i, sign, j) for j in range(d) for sign in (1, -1) for i in range(5)]
    edges = []
    for j in range(d):
        for i in range(5):
            next_j = (j + 1) % d if i == 4 else j
            edges.append(((i, 1, j), ((i + 1) % 5, 1, next_j), scale * (1 + eps * A[i])))
            edges.append(((i, -1, j), ((i + 1) % 5, -1, j), scale * (1 - eps * A[i])))
            edges.append(((i, 1, j), (i, -1, j), scale))
    return Graph(vertices, edges)


def fidelity():
    # Proposition 7.1(1): every equilateral trivalent cover has model balls.
    g = family(1, F(0), F(1))
    equilateral = 0
    for scale in (F(1, 3), F(1), F(7, 2)):
        h = family(1, F(0), scale)
        for v in h.vertices:
            for n in range(29):
                radius = scale * F(n, 8)
                assert h.ball(v, radius) == model(scale, radius)
                equilateral += 1
    # Example 7.2: x <= y <= z, x+y+z=3, z <= 3/2.
    theta = 0
    for ix in range(1, 25):
        for iy in range(ix, 25):
            x, y = F(ix, 8), F(iy, 8)
            z = 3 - x - y
            if not y <= z <= F(3, 2):
                continue
            h = Graph((0, 1), [(0, 1, x), (0, 1, y), (0, 1, z)])
            assert h.ball(0, z) == 6 * z - 3 == model(F(1), z)
            # Both functions are affine between these breakpoints on [0,z].
            points = sorted({F(0), x, y, F(1), z})
            for r in points:
                if 0 <= r <= z:
                    assert h.ball(0, r) >= model(F(1), r)
            theta += 1
    print(f'FIDELITY: equilateral vertex-radius-scale checks={equilateral}; theta metrics={theta}')
    return g


def affine_add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def affine_sum(xs):
    answer = (F(0), F(0))
    for x in xs:
        answer = affine_add(answer, x)
    return answer


def symbolic_local_identity(g, v, completed_levels, radius):
    """Exact coefficients of 1 and eps; cutoffs certified for all 0<eps<1/60."""
    def edge_affine(v, w):
        i, sign, _ = v
        k, other_sign, _ = w
        if sign != other_sign:
            return (F(1), F(0))
        edge_i = i if (i + 1) % 5 == k else k
        return (F(1), F(sign * A[edge_i]))

    volume = (F(0), F(0))
    stack = [(v, None, 0, (F(0), F(0)))]
    while stack:
        here, incoming, depth, distance = stack.pop()
        for e, w, _ in g.adj[here]:
            if e == incoming:
                continue
            length = edge_affine(here, w)
            end = affine_add(distance, length)
            if depth < completed_levels:
                # Endpoint distance < radius throughout the open eps interval.
                assert max(end[0], end[0] + end[1] * F(1, 60)) <= radius
                volume = affine_add(volume, length)
                stack.append((w, e, depth + 1, end))
            else:
                assert min(end[0], end[0] + end[1] * F(1, 60)) >= radius
                assert max(distance[0], distance[0] + distance[1] * F(1, 60)) <= radius
                volume = affine_add(volume, (radius - distance[0], -distance[1]))
    return volume


def settlement():
    assert tuple(A[i] + A[(i - 1) % 5] for i in range(5)) == T
    assert tuple(T[i] + T[(i - 1) % 5] + T[(i + 1) % 5] for i in range(5)) == U
    assert all(t * u < 0 for t, u in zip(T, U))
    assert model(F(1), F(3, 2)) == 6 and model(F(1), F(5, 2)) == 15
    assert F(7, 6) < F(3, 2) < 2 * F(5, 6)
    assert 2 * F(7, 6) < F(5, 2) == 3 * F(5, 6)
    symbolic_vertices = 0
    checked_vertices = 0
    ds = tuple(range(1, 13)) + (17, 31)
    for d in ds:
        base = family(d, F(0), F(1))
        assert base.connected()
        assert len(base.vertices) == 10 * d and len(base.edges) == 15 * d
        assert len(base.edges) - len(base.vertices) + 1 == 5 * d + 1
        assert all(v != w for v, w, _ in base.edges)
        assert len({frozenset((v, w)) for v, w, _ in base.edges}) == len(base.edges)
        assert all(len(base.adj[v]) == 3 for v in base.vertices)
        for v in base.vertices:
            i, sign, _ = v
            assert symbolic_local_identity(base, v, 1, F(3, 2)) == (F(6), F(-sign * T[i]))
            assert symbolic_local_identity(base, v, 2, F(5, 2)) == (F(15), F(-sign * U[i]))
            symbolic_vertices += 1
        for eps in (F(1, 1000), F(1, 120), F(59, 3600)):
            for scale in (F(1, 3), F(1), F(7, 2)):
                g = family(d, eps, scale)
                assert sum((e[2] for e in g.edges), F(0)) == 15 * d * scale
                assert all(F(5, 6) * scale < e[2] < F(7, 6) * scale for e in g.edges)
                for v in g.vertices:
                    i, sign, _ = v
                    s = sum((edge[2] for edge in g.adj[v]), F(0))
                    c = 2 * s + sum((sum((e[2] for e in g.adj[w]), F(0)) for _, w, _ in g.adj[v]), F(0))
                    r1, r2 = scale * F(3, 2), scale * F(5, 2)
                    b1, b2 = g.ball(v, r1), g.ball(v, r2)
                    assert s == scale * (3 + sign * eps * T[i])
                    assert b1 == 6 * r1 - s == scale * (6 - sign * eps * T[i])
                    assert b2 == 12 * r2 - c == scale * (15 - sign * eps * U[i])
                    assert b1 < model(scale, r1) or b2 < model(scale, r2)
                    checked_vertices += 1
    print('LOCAL: a=' + str(A) + '; t=' + str(T) + '; u=' + str(U))
    print('CUTOFFS: edge lengths in (5/6,7/6); radii 3/2 and 5/2; t_i*u_i<0 for every i')
    print(f'SYMBOLIC: exact affine identities on entire eps interval; vertices={symbolic_vertices}')
    print(f'FAMILY: d={ds}; eps=(1/1000,1/120,59/3600); scales=(1/3,1,7/2); vertices={checked_vertices}')
    print('BOUNDARY: eps=0 gives equality; positive eps excludes every lifted vertex at one of two radii')


if __name__ == '__main__':
    fidelity()
    settlement()
    print('ALL_CHECKS_OK')
