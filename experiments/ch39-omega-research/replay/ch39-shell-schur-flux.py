"""Flux of the exact seven-site Schur complement on the two clock torus."""
from functools import lru_cache
from pathlib import Path
import runpy

from flint import acb, acb_mat, arb

geometry = runpy.run_path(str(Path(__file__).with_name("ch39-shell-chiral-flux.py")))
shell = geometry["s"]
wall = geometry["wall"]
inside = tuple(i for i in range(19) if i not in wall)
modes = [(0, 0), (1, 0), (0, 1)]
harmonics = [(1, 0), (0, 1), (1, -1)]
pi = arb.pi()


def block(matrix, rows, columns):
    return acb_mat([[matrix[i, j] for j in columns] for i in rows])


def schur_data(u):
    matrix = shell["matrix"](u, {})
    interior = block(matrix, inside, inside)
    # The clock edges all point from the interior into the wall in the graph ordering.
    # Sylvester's criterion certifies that the conditional Gaussian exists.
    for n in range(1, len(inside) + 1):
        assert block(interior, range(n), range(n)).det().real > 0
    mixed = block(matrix, wall, inside)
    pieces = [mixed, acb_mat(len(wall), len(inside)), acb_mat(len(wall), len(inside))]
    for i, j, distance, _ in shell["edges"]:
        if distance in (724, 804):
            assert i in inside and j in wall
            row, col = wall.index(j), inside.index(i)
            pieces[1 + (distance == 804)][row, col] = mixed[row, col]
            mixed[row, col] = 0
    inverse = interior.inv()
    coefficients = {(0, 0): -block(matrix, wall, wall)}
    for h, left in zip(modes, pieces):
        for k, right in zip(modes, pieces):
            key = k[0] - h[0], k[1] - h[1]
            term = left * inverse * right.conjugate().transpose()
            coefficients[key] = coefficients.get(key, acb_mat(len(wall), len(wall))) + term
    return coefficients


def bond(coefficients, i, j):
    return {k: v[wall.index(i), wall.index(j)] for k, v in coefficients.items()}


def evaluate(poly, p, q):
    return sum(z * acb(0, a*p + b*q).exp() for (a, b), z in poly.items())


def imag_terms(poly, factor):
    coeff = {key: z * factor for key, z in poly.items()}
    terms = [(coeff[0, 0].imag, arb(0))]
    for a, b in harmonics:
        pos, neg = coeff[a, b], coeff[-a, -b]
        terms.append(((pos + neg).imag, (pos - neg).real))
    return terms


def area_enclosure(tests, depth):
    # Each box is in turns. The p-q interval has twice the radius of p or q.
    @lru_cache(maxsize=65536)
    def trig(level, center, width):
        scale = 1 << level
        center %= 2 * scale
        angle = pi * (center + arb(0, width)) / scale
        return angle.cos(), angle.sin()

    count_in = count_out = 0
    todo = [(0, 0, 0)]
    while todo:
        level, i, j = todo.pop()
        terms = (trig(level, 2*i+1, 1), trig(level, 2*j+1, 1),
                 trig(level, 2*(i-j), 2))
        values = [t[0][0] + sum(a*c + b*s for (a, b), (c, s) in zip(t[1:], terms))
                  for t in tests]
        area_units = 4**(depth-level)
        if all(v > 0 for v in values):
            count_in += area_units
        elif any(v < 0 for v in values):
            count_out += area_units
        elif level < depth:
            todo.extend((level+1, 2*i+di, 2*j+dj) for di in (0, 1) for dj in (0, 1))
    return arb(count_in) / 4**depth, 1 - arb(count_out) / 4**depth


def direct_degree(n, coefficients):
    h = acb.zeta_zero(n).imag
    p, q = shell["a"] * h, shell["b"] * h
    matrix = shell["matrix"](0, shell["phases"](h))
    effective = (block(matrix, wall, wall) - block(matrix, wall, inside)
                 * block(matrix, inside, inside).solve(block(matrix, inside, wall)))
    total = arb(0)
    for face in geometry["big"]:
        loop = acb(1)
        for i, j in zip(face, face[1:] + face[:1]):
            z = -effective[wall.index(i), wall.index(j)]
            assert (z - evaluate(bond(coefficients, i, j), p, q)).contains(0)
            loop *= z
        assert not loop.imag.contains(0) or loop.real > 0
        total += loop.arg()
    return int((total / (2*pi)).unique_fmpz())


if __name__ == "__main__":
    coefficients = schur_data(0)
    x = bond(coefficients, 7, 6)
    base = bond(coefficients, 11, 12)
    assert all(z == 0 for key, z in base.items() if key != (0, 0))
    y = base[0, 0]
    assert 0 < y.arg() < pi/3
    for pair in ((6, 14), (7, 12)):
        assert abs(evaluate(bond(coefficients, *pair), 0, 0)) > 0
    for n, degree in [(1, 0), (5, 3), (6, 3), (35, 0)]:
        assert direct_degree(n, coefficients) == degree
    low, high = area_enclosure([imag_terms(x, 1), imag_terms(x, -y)], 13)
    print("base bond phase", y.arg().str(15))
    print("Haar measure and asymptotic zero density", low.str(15), high.str(15))
    assert arb(".16") < low < high < arb(".18")
