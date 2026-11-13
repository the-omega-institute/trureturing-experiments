#!/usr/bin/env python3
"""Exact checks for the ternary (q,n,N)=(3,4,7) refutation candidate.

All evidence is integer or Fraction arithmetic.  The finite clique certificate
checks the structural lemma after normalizing a distance-four pair.
"""

from fractions import Fraction
from itertools import combinations, product


Word = tuple[int, ...]
WORDS: tuple[Word, ...] = tuple(product(range(3), repeat=4))


def distance(x: Word, y: Word) -> int:
    return sum(a != b for a, b in zip(x, y))


def volumes() -> tuple[int, ...]:
    return tuple(sum(2**j * __import__("math").comb(4, j) for j in range(t + 1)) for t in range(5))


def ball(center: Word, radius: int) -> set[Word]:
    return {x for x in WORDS if distance(center, x) <= radius}


def discrepancy(code: tuple[Word, ...]) -> Fraction:
    n = len(code)
    qn = len(WORDS)
    result = Fraction(0)
    vs = volumes()
    for t, v in enumerate(vs):
        for center in WORDS:
            m = sum(word in ball(center, t) for word in code)
            result += (Fraction(m, n) - Fraction(v, qn)) ** 2
    return result


def squared_count_sums(code: tuple[Word, ...]) -> tuple[int, ...]:
    return tuple(sum(sum(word in ball(center, t) for word in code) ** 2 for center in WORDS) for t in range(5))


def pair_intersections(x: Word, y: Word, t: int) -> int:
    return len(ball(x, t) & ball(y, t))


def max_clique(vertices: tuple[Word, ...]) -> tuple[Word, ...]:
    """Exact maximum clique in the distance-at-least-three graph."""
    neighbors = {
        v: {w for w in vertices if w != v and distance(v, w) >= 3}
        for v in vertices
    }
    best: list[Word] = []

    def expand(clique: list[Word], candidates: set[Word]) -> None:
        nonlocal best
        if len(clique) + len(candidates) <= len(best):
            return
        if not candidates:
            if len(clique) > len(best):
                best = clique[:]
            return
        # A deterministic branch order makes the certificate reproducible.
        while candidates:
            if len(clique) + len(candidates) <= len(best):
                return
            v = min(candidates)
            candidates.remove(v)
            expand(clique + [v], candidates & neighbors[v])

    expand([], set(vertices))
    return tuple(best)


def main() -> None:
    vs = volumes()
    assert vs == (1, 9, 33, 65, 81)
    # A perfect radius-e code would have 7 V_e = 81; this fails for every e.
    assert all(7 * v != 81 for v in vs)

    one_error_code = tuple(
        tuple(map(int, word))
        for word in ("0000", "0112", "0221", "1011", "1120", "1202", "2022")
    )
    distances = [distance(x, y) for x, y in combinations(one_error_code, 2)]
    assert len(distances) == 21 and set(distances) == {3}
    assert min(distances) == 3
    # Length four makes minimum distance at least five impossible for two words,
    # so the maximum correctable error count for seven-word codes is e=1.
    assert 2 * 2 + 1 > 4

    zero = (0, 0, 0, 0)
    one = (1, 1, 1, 1)
    normalized = tuple(w for w in WORDS if distance(w, zero) >= 3 and distance(w, one) >= 3)
    clique = max_clique(normalized)
    assert len(normalized) == 21 and len(clique) == 4
    # Therefore any distance-three code containing a distance-four pair has at
    # most 2+4=6 words.  Every seven-word e=1 code has all pair distances 3.

    candidate = tuple(
        tuple(map(int, word))
        for word in ("0000", "0012", "0221", "1111", "1202", "2122", "2210")
    )
    histogram = tuple(sum(distance(x, y) == d for x, y in combinations(candidate, 2)) for d in range(1, 5))
    assert histogram == (0, 1, 18, 2)
    assert discrepancy(one_error_code) == Fraction(11158, 3969)
    assert discrepancy(candidate) == Fraction(10348, 3969)
    assert discrepancy(one_error_code) - discrepancy(candidate) == Fraction(10, 49)
    assert squared_count_sums(one_error_code) == (7, 63, 735, 2597, 3969)
    assert squared_count_sums(candidate) == (7, 67, 721, 2597, 3969)

    expected_intersections = (
        (1, 0, 0, 0, 0),
        (9, 3, 2, 0, 0),
        (33, 21, 17, 12, 6),
        (65, 57, 53, 51, 50),
        (81, 81, 81, 81, 81),
    )
    for t in range(5):
        assert tuple(pair_intersections(zero, w, t) for w in (zero, (1, 0, 0, 0), (1, 1, 0, 0), (1, 1, 1, 0), one)) == expected_intersections[t]

    print("volumes", vs)
    print("one-error distance histogram", {d: distances.count(d) for d in sorted(set(distances))})
    print("normalized distance-four-pair candidate count", len(normalized), "maximum additional clique", len(clique))
    print("candidate distance histogram d=1..4", histogram)
    print("D(one-error code)", discrepancy(one_error_code), "D(candidate)", discrepancy(candidate), "gap", discrepancy(one_error_code) - discrepancy(candidate))
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
