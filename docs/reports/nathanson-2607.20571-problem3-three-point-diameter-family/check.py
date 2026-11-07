#!/usr/bin/env python3
"""Exact finite checks for Nathanson, arXiv:2607.20571v1, Problem 3.

All intervals and cubes are inclusive. No affine-rank or gcd counting formula
is used by the exhaustive enumeration. Translation removes duplicate triples.
"""
from itertools import combinations, product
from math import comb, gcd

H_MAX = 6


def sumsets_through(points, h_max):
    """Definition of hA: recursively add one element, with repetition allowed."""
    n = len(points[0])
    result = {(0,) * n}
    for _ in range(h_max):
        result = {tuple(x + y for x, y in zip(s, a))
                  for s in result for a in points}
        yield result


def freiman_fidelity():
    """Finite instances of source Theorem 1 and its base-g Freiman map."""
    cases = 0
    for n, d in ((2, 2), (3, 1)):
        base = H_MAX * d + 1  # Strictly greater than every sum coordinate.
        weights = [base ** i for i in range(n)]
        encode = lambda a: sum(x * w for x, w in zip(a, weights))
        points = list(product(range(d + 1), repeat=n))
        for triple in combinations(points, 3):
            integer_points = [(encode(a),) for a in triple]
            for lattice_sums, integer_sums in zip(
                    sumsets_through(triple, H_MAX),
                    sumsets_through(integer_points, H_MAX)):
                images = {encode(a) for a in lattice_sums}
                assert images == {a[0] for a in integer_sums}
                assert len(lattice_sums) == len(images)
                cases += 1
    print(f"FIDELITY: source Theorem 1 Freiman identity and injectivity; {cases} cases OK")


def triangular(r):
    return (r + 1) * (r + 2) // 2 if r >= 0 else 0


def published_spectrum(h):
    """Nathanson's published Theorem 9 (2411.02365v5), used as prior art."""
    return {triangular(h) - triangular(h - b) for b in range(2, h + 2)}


def published_lemma_check():
    cases = 0
    for b in range(2, 13):
        for a in range(1, b):
            if gcd(a, b) != 1:
                continue
            for h in range(1, H_MAX + 1):
                sums = {i * a + j * b for i in range(h + 1)
                        for j in range(h + 1 - i)}
                assert len(sums) == triangular(h) - triangular(h - b)
                cases += 1
    print(f"PRIOR_LEMMA: published three-point count; {cases} cases OK")


def exhaustive_cube_spectra(n, d_max):
    """Enumerate every translation class of triples in [0,d_max]^n.

    A triple is translated by its coordinatewise minimum to a unique triple
    whose coordinatewise minima are all zero. It fits in [0,D]^n exactly when
    its maximum coordinate after that translation is at most D. We enumerate
    all unordered triples in the largest cube, select these representatives,
    and calculate their actual sums. This exhausts every smaller cube too.

    Base H_MAX*d_max+1 is injective on all sum cubes. Subtracting h times the
    first encoded point leaves the same cardinality, allowing integer sets
    {i*u+j*v : i,j >= 0, i+j <= h}. No rank classification is assumed.
    """
    points = list(product(range(d_max + 1), repeat=n))
    base = H_MAX * d_max + 1
    weights = [base ** i for i in range(n)]
    encoded = {a: sum(x * w for x, w in zip(a, weights)) for a in points}
    coefficients = {h: [(i, j) for i in range(h + 1)
                        for j in range(h + 1 - i)]
                    for h in range(1, H_MAX + 1)}
    exact = [[set() for _ in range(H_MAX + 1)] for _ in range(d_max + 1)]
    representatives = [0] * (d_max + 1)
    raw_seen = 0
    for a, b, c in combinations(points, 3):
        raw_seen += 1
        if any(min(a[i], b[i], c[i]) != 0 for i in range(n)):
            continue
        d = max(max(a), max(b), max(c))
        representatives[d] += 1
        u, v = encoded[b] - encoded[a], encoded[c] - encoded[a]
        for h, pairs in coefficients.items():
            cardinality = len({i * u + j * v for i, j in pairs})
            exact[d][h].add(cardinality)
    assert raw_seen == comb(len(points), 3)
    # Independent inclusion-exclusion count of translation representatives.
    expected_reps = sum((-1) ** j * comb(n, j)
                        * comb(d_max ** j * (d_max + 1) ** (n - j), 3)
                        for j in range(n + 1))
    assert sum(representatives) == expected_reps
    accumulated = [set() for _ in range(H_MAX + 1)]
    spectra = {}
    for d in range(d_max + 1):
        for h in range(1, H_MAX + 1):
            accumulated[h].update(exact[d][h])
            spectra[d, h] = frozenset(accumulated[h])
            assert accumulated[h] <= published_spectrum(h)
        compact = "; ".join(f"h={h}:{sorted(accumulated[h])}"
                            for h in range(1, H_MAX + 1))
        print(f"CUBE n={n} D={d} all_triples={comb((d+1)**n, 3)} "
              f"translation_reps={sum(representatives[:d+1])} | {compact}")
    minima = {}
    for h in range(1, H_MAX + 1):
        minima[h] = next(d for d in range(1, d_max + 1)
                         if spectra[d, h] == published_spectrum(h))
    return spectra, minima


def main():
    freiman_fidelity()
    published_lemma_check()
    one_spectra, one_minima = exhaustive_cube_spectra(1, H_MAX + 1)
    for h in range(1, H_MAX + 1):
        assert one_minima[h] == (2 if h == 1 else h + 1)
    for n in (2, 3):
        spectra, minima = exhaustive_cube_spectra(n, H_MAX)
        for h in range(1, H_MAX + 1):
            assert minima[h] == (1 if h == 1 else h)
            if h >= 2:
                assert triangular(h) - 1 not in spectra[h - 1, h]
                # Equivalent to N_n(h,3) > N(h,3)^(1/n), with positive bases.
                assert minima[h] ** n > one_minima[h]
                print(f"REFUTATION h={h} n={n}: N={one_minima[h]} "
                      f"N_n={minima[h]}; {minima[h]}^{n}="
                      f"{minima[h] ** n}>{one_minima[h]}")
    print("MINIMA h=1..6: N(h,3)=" + str([one_minima[h] for h in range(1, 7)]))
    print("MINIMA h=1..6: N_2(h,3)=N_3(h,3)=[1, 2, 3, 4, 5, 6]")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
