"""Cellular chiral operator for three exact covers of the permutation shell."""
from itertools import combinations, permutations

from sympy import eye, ones, symbols, zeros


def shell_cells():
    vertices = list(permutations((0, 1, 3, 4)))
    # A permutohedron edge exchanges consecutive values in their sorted order.
    edges = []
    for i, j in combinations(range(24), 2):
        left, right = vertices[i], vertices[j]
        moved = tuple(sorted(x for x, y in zip(left, right) if x != y))
        if moved in ((0, 1), (1, 3), (3, 4)):
            edges.append((i, j))
    faces = [{j for j, p in enumerate(vertices) if p[i] == value}
             for value in (0, 4) for i in range(4)]
    faces += [{v for v, p in enumerate(vertices) if p[i]+p[j] == 1}
              for i, j in combinations(range(4), 2)]
    assert len(edges) == 36 and sorted(map(len, faces)) == [4]*6+[6]*8
    return vertices, edges, faces


def face_boundary(edges, faces):
    boundary = zeros(len(edges), len(faces))
    edge_ids = {e: i for i, e in enumerate(edges)}
    for f, face in enumerate(faces):
        start = min(face)
        previous, current = -1, start
        visited = set()
        while current not in visited:
            visited.add(current)
            neighbors = [v for v in face if tuple(sorted((current, v))) in edge_ids]
            assert len(neighbors) == 2
            dest = min(v for v in neighbors if v != previous)
            edge = tuple(sorted((current, dest)))
            boundary[edge_ids[edge], f] = 1 if dest > current else -1
            previous, current = current, dest
        assert current == start and visited == face
    return boundary


def exact_covers(faces):
    # Pick the least uncovered vertex and test just the facets through it.
    def extend(used, selected):
        if len(used) == 24:
            yield tuple(sorted(selected))
            return
        v = next(i for i in range(24) if i not in used)
        for j, face in enumerate(faces):
            if v in face and used.isdisjoint(face):
                yield from extend(used | face, selected+[j])
    return sorted(extend(set(), []))


def main():
    vertices, edges, faces = shell_cells()
    full = face_boundary(edges, faces)
    top, = full.nullspace()
    assert all(abs(v) == 1 for v in top)
    covers = exact_covers(faces)
    assert covers == [tuple(range(4)), tuple(range(4, 8)), tuple(range(8, 14))]
    shortest = {i for i,(v,w) in enumerate(edges)
                if sum((a-b)**2 for a,b in zip(vertices[v],vertices[w])) == 2}
    square_edges = {i for i,e in enumerate(edges) if any(set(e) <= f for f in faces[8:])}
    assert shortest == square_edges and len(shortest) == 24
    for cover in covers:
        labels = {v: j for j, f in enumerate(cover) for v in faces[f]}
        eids = [i for i, (v,w) in enumerate(edges) if labels[v] != labels[w]]
        fids = [i for i in range(14) if i not in cover]
        d2 = full.extract(eids, fids)
        d3 = top.extract(fids, [0])
        assert d2*d3 == zeros(len(eids), 1)
        # Q maps even forms to odd forms. Its adjoint is the other chiral block.
        Q = d2.col_join(d3.T)
        assert Q.rank() == Q.cols
        assert Q.rows-Q.cols == len(cover)-1
        transfer = zeros(len(eids), len(cover))
        for j, ei in enumerate(eids):
            v,w = edges[ei]
            transfer[j, labels[v]], transfer[j, labels[w]] = -1, 1
        assert d2.T*transfer == zeros(d2.cols, len(cover))
        assert transfer.rank() == len(cover)-1
        costs = [sum((a-b)**2 for a,b in zip(vertices[edges[i][0]], vertices[edges[i][1]]))//2
                 for i in eids]
        x = symbols("x")
        characteristic = (Q.T*Q).charpoly(x).as_expr()
        if len(cover) == 4:
            assert costs == [1]*12
            assert transfer.T*transfer == 8*eye(4)-2*ones(4)
            assert d2.rank() == 9 and d3.rank() == 1
            expected = (x-10)*(x-5)*(x-4)**3*(x-2)**2*(x-1)**3
        else:
            assert costs == [4]*12
            expected = (x-8)*(x-6)*(x-4)**3*(x-2)**3
        assert (characteristic-expected).expand() == 0
        print("cover", cover, "odd/even", Q.shape, "positive squared masses",
              characteristic.factor())


if __name__ == "__main__":
    main()
