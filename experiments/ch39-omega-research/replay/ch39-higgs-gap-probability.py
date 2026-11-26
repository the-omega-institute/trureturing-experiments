"""Exact constants for the finite Higgs and chiral spectral probability bounds."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from math import factorial


def kuhn_counts(side):
    # A triangle from its minimal vertex has two disjoint nonempty step masks.
    masks = range(1, 16)
    incidence = Counter()
    face_types = 0
    for first, second in product(masks, repeat=2):
        if first & second:
            continue
        incidence.update([first, second, first | second])
        face_types += 1
    assert face_types == 50 and len(incidence) == 15 and max(incidence.values()) == 14
    sites = side**4
    return sites, sites * len(incidence), sites * face_types, max(incidence.values())


def angular_bounds(degree):
    truncated_error = Fraction(2**23 * 16**97, factorial(97))
    # Use a short rational majorant before taking powers.
    assert truncated_error < Fraction(1, 10000)
    torus_weights = (1, 0, 0, -1)
    vandermonde_powers = [sum(a*b for a, b in zip(torus_weights, perm))
                          for perm in permutations(range(4))]
    assert max(vandermonde_powers) - min(vandermonde_powers) == 6
    fourier_degree = 96*degree + 6
    density = (2*fourier_degree + 1)*Fraction(10001, 9999)**degree
    assert density < 4096 and 4**2 * density < 256**2
    widths = []
    for exterior_power in (1, 2):
        charges = Counter(map(sum, combinations(torus_weights, exterior_power)))
        assert set(charges) == {-1, 0, 1} and charges[1] == charges[-1]
        # Four directed path occurrences, four spin components, both Laurent ends.
        widths.append(2*4*4*charges[1])
    assert widths == [32, 64]
    return max(widths)


def certify():
    sites, edges, faces, degree = kuhn_counts(5)
    width = angular_bounds(degree)
    mass = Fraction(139, 40)
    free_gap = Fraction(1, 64)
    assert (25 - 4*(mass-free_gap))**2 < 125 and Fraction(7, 2)-mass > free_gap
    assert abs(5-mass) + 10 < 16
    assert 2*(4+4)*2 == 32  # path length, Clifford plus Wilson, exterior square
    cotree = 4*sites - (sites-1)
    angles = 3*cotree
    loss = 12
    dim = 24*sites
    gap_bits = 10*dim + 2*width*(8*angles+loss)
    tube_bits = gap_bits + 8
    kinetic_bits = 2*tube_bits + 27
    assert (sites, edges, faces, cotree, angles) == (625, 9375, 31250, 1876, 5628)
    assert (gap_bits, tube_bits, kinetic_bits) == (5914608, 5914616, 11829259)
    # b < 2^-14 and 2D < 2^-14, retaining exponents symbolically.
    assert edges < 2**14 and 5 < 2**3
    assert 3*2**25 > 14 + 23*degree + 10*(tube_bits+3) + 14
    assert (140*faces).bit_length() + 14 < tube_bits
    assert Fraction(2, 2**loss) + Fraction(2, 2**14) < Fraction(1, 1000)
    assert 2*32 < 2**7
    # A larger gap at the previously certified kinetic coefficient.
    assert free_gap - Fraction(32, 2**13) - Fraction(64, 2**52) > Fraction(1, 100)
    assert 1+4*2**13 < 2**16
    positive_sector_bits = (23*degree+32*16)*cotree
    assert positive_sector_bits == 1564584
    assert 3*2**9 > 14+23*degree+10*(52+3)+14
    assert Fraction(140*faces, 2**52) < Fraction(1, 2**14)
    print(f"5^4: reduced gap >= 2^-{gap_bits}, soft kappa = 2^{kinetic_bits}.")
    print(f"At kappa = 2^115, P(g > 1/100 and section defined) > 2^-{positive_sector_bits+1}.")


if __name__ == "__main__":
    certify()
