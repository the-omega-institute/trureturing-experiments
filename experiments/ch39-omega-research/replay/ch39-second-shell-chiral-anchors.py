"""比较第二个十九点整数壳中两组三点近锚的极值三角与有向钟。"""

from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from runpy import run_path
import sys


def rotate(point):
    x, y, z = point
    return y, z, x


def reflect(point):
    x, y, z = point
    return y, x, z


def rotate_face(face):
    return tuple(sorted(rotate(point) for point in face))


def count_faces(faces, clocks):
    for clock in clocks:
        assert {(rotate(p), rotate(q)) for p, q in clock} == clock
    counts = Counter()
    orbit_counts = Counter()
    for face in faces:
        traversal = tuple(zip(face, face[1:] + face[:1]))
        freq = tuple(sum(int(edge in clock) - int(edge[::-1] in clock) for edge in traversal)
                     for clock in clocks)
        label = max(freq, tuple(-x for x in freq))
        orbit = {face, rotate_face(face), rotate_face(rotate_face(face))}
        assert orbit <= faces
        if label != (0, 0):
            assert len(orbit) == 3
        counts[label] += 1
        if min(orbit) == face:
            orbit_counts[label] += 1
    assert set(counts) == {(0, 0), (1, 0), (0, 1), (1, -1)}
    active = sum(counts.values()) - counts[0, 0]
    k = sum(orbit_counts.values()) - orbit_counts[0, 0]
    assert active == 3 * k
    return counts, k


def extreme_clocks(bases, anchors, distance):
    clocks = [set(), set()]
    for p, q in bases:
        candidates = {a: tuple(sorted((distance(p, a), distance(q, a)))) for a in anchors}
        assert sorted(candidates.values()) == [(216, 296), (536, 776), (776, 936)]
        apex, = [a for a, sides in candidates.items()
                 if all(sides[0] >= x and sides[1] >= y for x, y in candidates.values())]
        for clock, length in zip(clocks, candidates[apex]):
            source, = [z for z in (p, q) if distance(z, apex) == length]
            clock.add((source, apex))
    return clocks


def verify(root):
    lattice = run_path(str(root / "ch39-nineteen-lattice-shells.py"))
    points = sorted(lattice["shell"](1083))
    distance = lattice["square_distance"]
    center = (19, 19, 19)
    radii = {p: distance(p, center) for p in points if p != center}
    near = {p for p, r in radii.items() if r == min(radii.values())}
    far = {p for p, r in radii.items() if r == max(radii.values())}
    first = min(near)
    anchors = {first, rotate(first), rotate(rotate(first))}
    opposite = {reflect(p) for p in anchors}
    assert len(anchors) == len(opposite) == 3 and anchors.isdisjoint(opposite)
    assert anchors | opposite == near
    far_pairs = list(combinations(sorted(far), 2))
    base = min(distance(*pair) for pair in far_pairs)
    bases = [pair for pair in far_pairs if distance(*pair) == base]
    assert base == 200 and len(bases) == 3 and {p for pair in bases for p in pair} == far
    assert sorted(distance(*pair) for pair in bases) == [200] * 3

    edges = {e for e in combinations(points, 2) if distance(*e) <= 936}
    degrees = Counter(p for edge in edges for p in edge)
    faces = {tri for tri in combinations(points, 3) if all(e in edges for e in combinations(tri, 2))}
    assert len(edges) == 150 and len(faces) == 676
    assert Counter(degrees.values()) == {13: 6, 16: 6, 18: 7}
    assert {p for p in points if degrees[p] == 13} == far
    assert all(degrees[p] == 18 for p in near)
    ratio = Fraction(776, base)
    assert ratio == Fraction(97, 25) and (724 / ratio).denominator == 97
    print("近锚的两个保向轨道", sorted(anchors), sorted(opposite))
    print("底边/短侧/长侧平方、图边数、三角数、短侧与底边的平方比", (base, 776, 936),
          len(edges), len(faces), ratio)

    all_clocks = [{(p, q) for p in far for q in near if distance(p, q) == length}
                  for length in (776, 936)]
    for clock, length in zip(all_clocks, (776, 936)):
        assert {tuple(sorted(e)) for e in clock} == {e for e in edges if distance(*e) == length}
    marked = []
    for group in (anchors, opposite):
        picked = extreme_clocks(bases, group, distance)
        marked.append(picked)
        by_endpoint = [{e for e in clock if e[1] in group} for clock in all_clocks]
        assert all(c <= larger for c, larger in zip(picked, by_endpoint))
        for label, clocks, sizes, expected, expected_k in (
            ("极值三角实际侧边", picked, (3, 3), (607, 33, 33, 3), 23),
            ("同长且终点在选定近锚", by_endpoint, (6, 3), (580, 60, 30, 6), 32),
            ("全图同长边", all_clocks, (12, 6), (496, 108, 48, 24), 60),
        ):
            assert tuple(map(len, clocks)) == sizes
            counts, k = count_faces(faces, clocks)
            classes = tuple(counts[f] for f in ((0, 0), (1, 0), (0, 1), (1, -1)))
            assert classes == expected and k == expected_k
            print("近锚代表、定钟条件、两类边数、四类面数、K、D",
                  min(group), label, sizes, classes, k, k * ratio)
    for first_clock, other_clock in zip(*marked):
        assert first_clock.isdisjoint(other_clock)
        assert {(reflect(p), reflect(q)) for p, q in first_clock} == other_clock


if __name__ == "__main__":
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
    verify(root)
