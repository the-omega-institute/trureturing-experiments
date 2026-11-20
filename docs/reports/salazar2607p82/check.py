"""Exact checks of Salazar's GRWS samples, atomic representation, and increment bounds."""

from fractions import Fraction
from itertools import product


def sample(p, numerator_offset, denominator_offset, n):
    return (p**n + numerator_offset) / (p**n + denominator_offset)


def check_fidelity(p, numerator_offset, denominator_offset, n):
    q = 1 / p
    weight_total_factor = denominator_offset - numerator_offset
    decay = -denominator_offset
    z = q**n
    value = sample(p, numerator_offset, denominator_offset, n)
    boundary = sample(p, numerator_offset, denominator_offset, 0)
    assert value == (1 + numerator_offset * z) / (1 + denominator_offset * z)

    # Equation (84): the m-th atom is at m log(p), so its kernel at n is 1 - z**m.
    cutoff = 5
    masses = [weight_total_factor * decay ** (m - 1) for m in range(1, cutoff + 1)]
    partial = sum(mass * (1 - z**m) for m, mass in enumerate(masses, start=1))
    mass_tail = weight_total_factor * decay**cutoff / (1 - decay)
    weighted_tail = weight_total_factor * z * (decay * z) ** cutoff / (1 - decay * z)
    assert all(mass >= 0 for mass in masses)
    assert masses[0] > 0
    assert boundary + partial <= value <= boundary + partial + mass_tail
    assert value == boundary + partial + mass_tail - weighted_tail


def check_increment(p, numerator_offset, denominator_offset, n):
    q = 1 / p
    increment = (
        sample(p, numerator_offset, denominator_offset, n + 1)
        - sample(p, numerator_offset, denominator_offset, n)
    )
    factor = (denominator_offset - numerator_offset) * (1 - q)
    denominator = (1 + denominator_offset * q**n) * (1 + denominator_offset * q ** (n + 1))
    assert denominator >= (1 + denominator_offset) ** 2 > 0
    assert increment == factor * q**n / denominator
    assert 0 < increment <= factor * q**n / (1 + denominator_offset) ** 2


def check_interval_comparison(p, numerator_offset, denominator_offset, cut, density_floor):
    q = 1 / p
    upper_constant = (denominator_offset - numerator_offset) * (1 - q) / (1 + denominator_offset) ** 2

    # Let r_v = exp(-v) and r_u = exp(-u). Then v-u >= 1-r_v/r_u.
    r_v = q + (1 - q) * cut
    r_u = (1 + r_v) / 2
    assert 0 < q < r_v < r_u < 1
    lower_constant = (1 - r_v / r_u) * (1 - r_u) * density_floor
    assert lower_constant > 0
    ratio = r_v / q
    n = 1
    while lower_constant * ratio**n <= upper_constant:
        n *= 2
    assert lower_constant * r_v**n > upper_constant * q**n
    check_increment(p, numerator_offset, denominator_offset, n)
    for position in range(9):
        w = r_v + (r_u - r_v) * Fraction(position, 8)
        assert w**n * (1 - w) * density_floor >= r_v**n * (1 - r_u) * density_floor
    return n


def main():
    bases = tuple(Fraction(n, d) for n, d in [(8, 7), (5, 4), (3, 2), (2, 1), (5, 1)])
    numerators = (Fraction(-7, 8), Fraction(-1, 2), Fraction(-1, 8))
    denominators = (Fraction(-3, 4), Fraction(-1, 4), Fraction(-1, 16), Fraction(0))
    triples = [
        (p, numerator_offset, denominator_offset)
        for p, numerator_offset, denominator_offset in product(bases, numerators, denominators)
        if -1 < numerator_offset < denominator_offset <= 0
    ]
    indices = range(13)
    cuts = (Fraction(1, 4), Fraction(1, 2), Fraction(7, 8))
    density_floors = (Fraction(1, 1000), Fraction(1, 3), Fraction(5))
    thresholds = []
    for p, numerator_offset, denominator_offset in triples:
        for n in indices:
            check_fidelity(p, numerator_offset, denominator_offset, n)
            check_increment(p, numerator_offset, denominator_offset, n)
        for cut, density_floor in product(cuts, density_floors):
            thresholds.append(
                check_interval_comparison(p, numerator_offset, denominator_offset, cut, density_floor)
            )
    for p, offset, n in product(bases, denominators, indices):
        assert sample(p, offset, offset, n) == 1
    boundary_triples = sum(denominator_offset == 0 for _, _, denominator_offset in triples)
    print(f"Equations (78), (80), (84): {len(triples) * len(indices)} exact sample checks")
    print(f"Admissible triples: {len(triples)}, with D = 0 in {boundary_triples}")
    print(f"Increment identity and geometric upper bound: {len(triples) * len(indices)} exact checks")
    print(f"Interval lower bound crosses the upper bound in {len(thresholds)} exact checks")
    print(f"Decisive powers of two: {min(thresholds)} through {max(thresholds)}")
    print(f"N = D constant boundary: {len(bases) * len(denominators) * len(indices)} exact checks")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
