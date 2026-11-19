#!/usr/bin/env python3
"""Exact replica contractions from computational-basis support and torus gluing."""

from fractions import Fraction
from itertools import combinations, product


def sites(n, dimension):
    return list(product(range(n), repeat=dimension))


def shifts(n, site):
    yield site
    for axis in range(len(site)):
        yield site[:axis] + ((site[axis] + 1) % n,) + site[axis + 1 :]


def ghz_contraction_exponent(n, q):
    """Enumerate the ket support and test the shifted bra support in Eq. (5)."""
    torus = sites(n, q - 1)
    index = {site: i for i, site in enumerate(torus)}
    bra_positions = [tuple(index[s] for s in shifts(n, site)) for site in torus]
    count = 0
    for assignment in range(1 << len(torus)):
        if all(
            all((assignment >> r & 1) == (assignment >> positions[0] & 1)
                for r in positions[1:])
            for positions in bra_positions
        ):
            count += 1
    assert count == 2
    # Two nonzero amplitudes 1/sqrt(2) per site give weight 2^(-|torus|).
    return len(torus) - 1


def gf_multiply(a, b, modulus):
    """Multiply binary polynomials modulo a monic irreducible polynomial."""
    degree = modulus.bit_length() - 1
    value = 0
    while b:
        if b & 1:
            value ^= a
        b >>= 1
        a <<= 1
        if a & (1 << degree):
            a ^= modulus
    return value


def support_words(d, parameter, modulus):
    words = [
        (a, a ^ s, s, s ^ gf_multiply(parameter, a, modulus))
        for a, s in product(range(d), repeat=2)
    ]
    assert len(set(words)) == d**2
    assert Fraction(len(words), d**2) == 1
    for i, j in combinations(range(4), 2):
        assert len({(word[i], word[j]) for word in words}) == d**2
    return words


def pack(word, degree):
    return sum(component << (degree * party) for party, component in enumerate(word))


def support_checks(words, degree):
    """Recover the dual binary code from the literal support words."""
    packed = [pack(word, degree) for word in words]
    checks = [
        mask for mask in range(1, 1 << (4 * degree))
        if all((mask & word).bit_count() % 2 == 0 for word in packed)
    ]
    # Use an independent basis of the dual code in the replica matrix.
    basis = list(row_basis(checks).values())
    assert len(basis) == 2 * degree
    return basis


def row_basis(rows):
    pivots = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = row
                break
            row ^= pivots[pivot]
    return pivots


def replica_nullity(n, degree, checks):
    """Rank the ket and bra constraints on all four labels at each torus site."""
    torus = sites(n, 3)
    index = {site: i for i, site in enumerate(torus)}
    register_mask = (1 << degree) - 1

    def rows():
        for i, site in enumerate(torus):
            bra_sites = tuple(index[s] for s in shifts(n, site))
            for check in checks:
                yield check << (4 * degree * i)
                yield sum(
                    ((check >> (degree * p)) & register_mask)
                    << (degree * (4 * bra_site + p))
                    for p, bra_site in enumerate(bra_sites)
                )

    rank = len(row_basis(rows()))
    return 4 * degree * len(torus) - rank


def main():
    # Section 6, Table 2. The displayed entropies are in units of log(2).
    for q, n in ((2, 2), (2, 3), (2, 4), (3, 2), (3, 3), (3, 4)):
        exponent = ghz_contraction_exponent(n, q)
        entropy = Fraction(exponent, (n - 1) * n ** (q - 2))
        expected = Fraction(1) if q == 2 else 1 + Fraction(1, n)
        assert entropy == expected
        print(f"Table 2 q={q} n={n}: S/log(2)={entropy}")

    # The irreducible binary polynomials are x²+x+1 and x³+x+1.
    for modulus in (0b111, 0b1011):
        degree = modulus.bit_length() - 1
        d = 1 << degree
        for parameter in range(2, d):
            words = support_words(d, parameter, modulus)
            checks = support_checks(words, degree)
            if d == 4 and parameter == 2:
                for a, s in product(range(4), repeat=2):
                    a0, a1 = a & 1, a >> 1
                    s0, s1 = s & 1, s >> 1
                    bits = (a0, a1, a0 ^ s0, a1 ^ s1,
                            s0, s1, s0 ^ a1, s1 ^ a0 ^ a1)
                    assert pack(words[d * a + s], degree) == pack(bits, 1)
            entropies = []
            for n, expected_nullity in (
                (d, 3 * d**2 // 2), (d + 1, (d + 1)**2 + d + 1)
            ):
                nullity = replica_nullity(n, degree, checks)
                assert nullity == degree * expected_nullity
                exponent = 2 * degree * n**3 - nullity
                entropy = Fraction(exponent, (n - 1) * n**2)
                entropies.append(entropy)
                print(f"d={d} lambda={parameter} n={n}: "
                      f"binary nullity={nullity}, Z=2^(-{exponent}), S/log(2)={entropy}")
            gap = entropies[1] - entropies[0]
            assert gap == Fraction(degree * (d - 3), 2 * (d**2 - 1)) > 0
            print(f"d={d} lambda={parameter}: gap/log(2)={gap}")

    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
