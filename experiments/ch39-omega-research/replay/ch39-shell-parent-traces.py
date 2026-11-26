"""Restriction of three spin trace forms, and the mixed gravitational coefficient."""
from fractions import Fraction
from itertools import product


# (complex dimension, central charge, SU4, SU3, weak left, weak right).
# The four nonabelian indices are normalized to one on a defining block.
mixed = (12, 1, 3, 4, 0, 0)
left = (2, 0, 0, 0, 1, 0)
right = (2, 0, 0, 0, 0, 1)
internal = (mixed, (3, -2, 0, 1, 0, 0), (1, -6, 0, 0, 0, 0),
            (15, 0, 8, 0, 0, 0), (8, 0, 0, 6, 0, 0))
parents = [(mixed, left, right), (mixed, right, right),
           internal + (left,) * 16 + (right,) * 14]


def levels(blocks):
    return tuple(sum(b[i] for b in blocks) for i in range(2, 6))


def central_metric(blocks):
    return sum(Fraction(b[0] * b[1] ** 2, 144) for b in blocks)


def residual_cost(rank, flux):
    residue = flux % rank
    return Fraction(residue * (rank - residue), rank)


def energy(blocks, flux, face):
    offset_a, offset_b = (1, 0) if face == "cap" else (-1, -1)
    w4, w3, _, _ = levels(blocks)
    return (w4 * residual_cost(4, flux - offset_a)
            + w3 * residual_cost(3, -flux - offset_b)
            + central_metric(blocks) * flux ** 2)


assert [sum(b[0] for b in blocks) for blocks in parents] == [16, 16, 99]
# Every determinant is a power of J^2, including the trivial determinant.
assert [sum(b[0] * b[1] for b in blocks) for blocks in parents] == [12, 12, 0]
rows = [levels(blocks) for blocks in parents]
assert rows == [(3, 4, 1, 1), (3, 4, 0, 2), (11, 11, 16, 14)]
assert tuple(map(sum, zip(*rows))) == (17, 19, 17, 17)
assert tuple(a - b for a, b in zip(rows[1], rows[0])) == (0, 0, -1, 1)
assert [central_metric(blocks) for blocks in parents] == [Fraction(1,12), Fraction(1,12), Fraction(5,12)]
# The bare virtual bulk representation is E*L minus conjugate(E)*R.
# Reversing the complex representation preserves its quadratic index.
bare_indices = (2*3-2*3, 2*4-2*4, 12, -12)
assert bare_indices == (0, 0, 12, -12)

for third_weight, side_min in [(Fraction(0), [-5]), (Fraction(1, 4), [-5, -1]),
                               (Fraction(1, 2), [-1])]:
    for face, expected in (("cap", [-3]), ("side", side_min)):
        energies = {k: energy(parents[0], k, face) + third_weight * energy(parents[2], k, face)
                    for k in range(-5, 6, 2)}
        assert [k for k, e in energies.items() if e == min(energies.values())] == expected
indices = set()
for sides in product((-1, -5), repeat=6):
    total = -6 + sum(sides)
    assert total % 4 == 0
    indices.add(-total // 4)
assert indices == set(range(3, 10))
print("spin defining real dimensions after padding: 198, 198, 198")
print("equal kinetic weights give SU4, SU3-family, left, right levels:", *map(sum, zip(*rows)))
print("bare mixed gravitational signed quadratic indices:", *bare_indices)
print("coexisting face choices at the critical ratio give net indices:", sorted(indices))
