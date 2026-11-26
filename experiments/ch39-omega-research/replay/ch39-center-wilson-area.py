"""Center-character selection and finite filling-area bounds on SU(5) Wilson loops."""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from json import loads
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def tetrahedron():
    edges = list(combinations(range(4), 2))
    faces = list(combinations(range(4), 3))
    counts = [[0]*5 for _ in range(16)]
    for phase in product(range(5), repeat=6):
        curls = [(phase[edges.index((i, j))]+phase[edges.index((j, k))]-phase[edges.index((i, k))]) % 5
                 for i, j, k in faces]
        flat = sum(1 << i for i, v in enumerate(curls) if v == 0)
        mask = flat
        while True:
            counts[mask][curls[0]] += 1
            if not mask:
                break
            mask = (mask-1) & flat
    for mask, bins in enumerate(counts):
        spans = bool(mask & 1) or mask == 14
        assert bins[0] > 0
        if spans:
            assert not any(bins[1:])
        else:
            assert len(set(bins)) == 1


def coordinate_projection():
    types = set()
    for order in permutations(range(4)):
        vertex = (0,)*4
        vertices = [vertex]
        for axis in order:
            vertex = tuple(x+int(i == axis) for i, x in enumerate(vertex))
            vertices.append(vertex)
        for face in combinations(vertices, 3):
            triangle = tuple(v[:2] for v in face)
            if len(set(triangle)) == 3:
                types.add(triangle)
    assert types == {((0, 0), (1, 0), (1, 1)), ((0, 0), (0, 1), (1, 1))}


def verify():
    tetrahedron()
    coordinate_projection()
    cells = loads((ROOT/'ch39-tick-separators.json').read_text())['prefix']
    def rate(height):
        strength = 1/(height**2+F(1, 4))
        b = (30+F(2, 15))*strength
        return 56*(2*b/(1+b))
    assert rate(F(cells[11][0]+1, 1 << 400)) > 1
    for j, upper in [(13, F('.951')), (46, F('.186')), (3200, F(1, 4096))]:
        value = rate(F(cells[j-1][0], 1 << 400))
        assert 0 < value < upper < 1
        print('First allowed zero, filling-area factor:', j, float(value))


if __name__ == '__main__':
    verify()
