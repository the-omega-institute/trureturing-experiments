"""Certify the exchange selection rule for arbitrary left and right shell currents."""
from functools import lru_cache
from itertools import combinations_with_replacement
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, arb_mat, ctx

root = Path(__file__).resolve().parent
model = run_path(str(root / 'ch39-shell-propagator-readings.py'))
native = model['native']
points = model['points']
orbits = model['groups']
indices = list(range(1, 31)) + [
    33, 38, 42, 47, 53, 60, 67, 76, 85, 96, 108, 121, 136, 153, 172, 193, 217, 244,
    274, 308, 347, 390, 438, 492, 554, 622, 699, 786, 884, 993, 1117, 1255, 1411,
    1586, 1783, 2004, 2253, 2533, 2847, 3200,
]


def rotate(point):
    return point[1:] + point[:1]


def reflect(point):
    return point[1], point[0], point[2]


@lru_cache(None)
def layers():
    diagonal, hopping = {}, {}
    size = len(points)
    for i, j in native['links']:
        x, y = points[i], points[j]
        length = sum((a - b)**2 for a, b in zip(x, y, strict=True))
        if length not in diagonal:
            diagonal[length] = acb_mat(size, size)
            hopping[length] = acb_mat(size, size)
        a, b = native['link'](i, j)
        for transform, expected in [(rotate, (a, b)), (reflect, (a, -b))]:
            edge = points.index(transform(x)), points.index(transform(y))
            assert edge in native['links'] and native['link'](*edge) == expected
        phase = acb(a, arb(3).sqrt() * b) / arb(a*a + 3*b*b).sqrt()
        diagonal[length][i, i] += 1
        hopping[length][i, j] = -phase
    assert len(diagonal) == 23
    return diagonal, hopping


@lru_cache(None)
def sector_basis(k):
    # Fourier orbit vectors fixed by reflection followed by complex conjugation.
    columns = []
    taken = set()
    omega = acb(-1, arb(3).sqrt()) / 2
    for orbit in orbits:
        if orbit[0] in taken:
            continue
        if len(orbit) == 1:
            if k == 0:
                v = acb_mat(19, 1)
                v[orbit[0], 0] = 1
                columns.append(v)
            continue
        fixed = [points[i] for i in orbit if reflect(points[i]) == points[i]]
        point = fixed[0] if fixed else points[orbit[0]]
        v, reflected = acb_mat(19, 1), acb_mat(19, 1)
        for step in range(3):
            weight = omega**((k * step) % 3) / arb(3).sqrt()
            i, j = points.index(point), points.index(reflect(point))
            v[i, 0], reflected[j, 0] = weight, weight.conjugate()
            taken.update((i, j))
            point = rotate(point)
        if fixed:
            columns.append(v)
        else:
            columns.extend([(v + reflected) / arb(2).sqrt(),
                            (v - reflected) * acb(0, 1) / arb(2).sqrt()])
    return acb_mat([[v[i, 0] for v in columns] for i in range(19)])


@lru_cache(None)
def sector_layers(k):
    basis = sector_basis(k)
    star = basis.conjugate().transpose()
    result = []
    for full in layers():
        reduced = {}
        for length, layer in full.items():
            block = star * layer * basis
            assert all(v.imag.contains(0) for v in block.entries())
            reduced[length] = arb_mat([[v.real for v in r] for r in block.tolist()])
        result.append(reduced)
    return result


def edge_coefficients(s, continuous):
    terms = {}
    for length in layers()[0]:
        if continuous:
            t = (s * model['density'](length)).sqrt()
            terms[length] = t / t.tanh(), t / t.sinh()
        else:
            terms[length] = arb(1), arb(1)
    return terms


def green(s, terms, diagonal, hopping, matrix_type):
    size = next(iter(diagonal.values())).ncols()
    pencil = matrix_type([[s if i == j else 0 for j in range(size)] for i in range(size)])
    for length, (a, b) in terms.items():
        pencil += a * diagonal[length] + b * hopping[length]
    return pencil.inv()


def exchange_row(s):
    terms = edge_coefficients(s, True)
    row = []
    for k in range(3):
        block = green(s, terms, *sector_layers(k), arb_mat)
        row.extend(block[i, j] for i, j in combinations_with_replacement(range(block.ncols()), 2))
    return row


def verify_span(rows, bound):
    matrix = arb_mat(rows)
    n = matrix.nrows()
    assert matrix.ncols() == n
    identity = arb_mat([[int(i == j) for j in range(n)] for i in range(n)])
    residual = matrix * matrix.mid().inv().mid() - identity
    norm = max(sum(x.abs_upper() for x in row).upper() for row in residual.tolist())
    assert norm < arb(bound)
    print(n, 'independent columns, maximum absolute row sum of residual', norm.str(7))


def main():
    with ctx.workprec(2048):
        heights = {i: acb.zeta_zero(i).imag for i in indices + [3199]}
        assert all(model['height'](i).contains(g) for i, g in heights.items())
    with ctx.workprec(4096):
        shifts = {i: (g / heights[1])**2 for i, g in heights.items()}
        assert [sector_basis(k).ncols() for k in range(3)] == [7, 6, 6]
        rows = [exchange_row(shifts[i]) for i in indices]
        verify_span(rows, '1e-428')
        verify_span([r + [arb(1)] for r in rows + [exchange_row(shifts[3199])]], '1e-428')
        for continuous in (False, True):
            responses = [green(shifts[i], edge_coefficients(shifts[i], continuous),
                               *layers(), acb_mat) for i in range(1, 20)]
            for orbit in orbits:
                if len(orbit) == 3:
                    samples = acb_mat([[g[orbit[0], j] for j in range(19)] for g in responses])
                    assert abs(samples.det()) > 0
            print('First nineteen detect all currents from each nonpolar site, arc interior', continuous)


if __name__ == '__main__':
    main()
