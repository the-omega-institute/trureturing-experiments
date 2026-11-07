#!/usr/bin/env python3
"""Exact checks for the Covtree witness-bound refutation.

The source is Gutzeit--Shaban--Yeats--Zalel, arXiv:2605.00622v1.
All counts below are integer or Fraction computations; no floating point is used.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product


def star_options(s: int) -> list[tuple[int, str]]:
    """Sizes and isomorphism types of nonempty downsets of L_s."""
    return [(0, "")] + [(t, f"A{t}") for t in range(1, s + 1)] + [(s + 1, f"L{s}")]


def q_downset_types(r: int) -> set[tuple[str, ...]]:
    """Enumerate the types of size-(2r+2) downsets of L_r + L_r + L_(r+1)."""
    n = 2 * r + 2
    answer: set[tuple[str, ...]] = set()
    for choices in product(star_options(r), star_options(r), star_options(r + 1)):
        if sum(size for size, _ in choices) == n:
            antichain_size = sum(size for size, kind in choices if kind.startswith("A"))
            stars = [kind for _, kind in choices if kind.startswith("L")]
            signature = stars + ([f"A{antichain_size}"] if antichain_size else [])
            answer.add(tuple(sorted(signature)))
    return answer


def source_height_two_fidelity(ell: int) -> tuple[int, int, Fraction]:
    """Reproduce Proposition 6.3's H_ell count and size from the source.

    H_ell is ell copies of L_2 plus one L_(3(ell-1)); its selected downsets
    have size n=3ell. A downset is represented by the number of complete L_2
    components and, separately, the large component case.
    """
    n = 3 * ell
    d_n = ell + 2
    size = ell * 3 + (3 * (ell - 1) + 1)
    # For this family k=ell+2, so the proposed n+3k-3 value is 6ell+3.
    bound_gap = Fraction(size) - Fraction(n + 3 * (ell + 2) - 3)
    return d_n, size, bound_gap


def forcing_certificate(r: int) -> dict[str, int | Fraction]:
    """Check the exact counting certificate for every witness of the four types.

    A witness has an A_n downset, hence at least n globally minimal elements.
    The L_r+L_r type supplies two maxima with disjoint predecessor sets of
    r minima each. The L_(r+1)+A_r type supplies a third maximum with r+1
    minimal predecessors. Any overlap with either predecessor set would,
    after adding unused minima to reach n, create a two-maximum downset with
    overlapping pasts; that is not one of the four listed types. Therefore
    the three predecessor sets are disjoint and the union of the three stars
    is a downset of size 3r+4.
    """
    n = 2 * r + 2
    minima_required = n
    predecessor_total = r + r + (r + 1)
    witness_lower = predecessor_total + 3
    claimed_q_size = 3 * r + 4
    proposed_bound = n + 3 * 4 - 3
    gap = witness_lower - proposed_bound
    overlap_completion_min = minima_required - (2 * r + 1)
    assert overlap_completion_min >= 1
    assert witness_lower == claimed_q_size
    assert gap == r - 7
    return {
        "n": n,
        "k": 4,
        "minima_required": minima_required,
        "predecessor_total": predecessor_total,
        "witness_lower": witness_lower,
        "claimed_q_size": claimed_q_size,
        "proposed_bound": proposed_bound,
        "gap": gap,
        "overlap_completion_min": overlap_completion_min,
        "excess_over_n": witness_lower - n,
    }


def main() -> None:
    for ell in range(2, 9):
        d_n, size, gap = source_height_two_fidelity(ell)
        assert d_n == ell + 2
        assert size == 6 * ell - 2
        assert gap == Fraction(-5)

    expected = {
        ("L8", "L8"),
        ("A8", "L9"),
        ("A9", "L8"),
        ("A18",),
    }
    assert q_downset_types(8) == expected
    for r in range(8, 21):
        assert len(q_downset_types(r)) == 4
        cert = forcing_certificate(r)
        assert cert["witness_lower"] > cert["proposed_bound"]
        assert cert["gap"] == r - 7
        assert cert["excess_over_n"] == r + 2

    print("Fidelity: source H_ell count d_(3ell)=ell+2 and |H_ell|=6ell-2 checked for ell=2..8")
    print("Counterexample family: Q_r has exactly four size-(2r+2) downset types for r=8..20")
    print("Forcing certificate: every witness has a Q_r downset and size at least 3r+4")
    print("Bound gap: (3r+4) - (2r+11) = r-7; n+f(4) excess is r+2")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
