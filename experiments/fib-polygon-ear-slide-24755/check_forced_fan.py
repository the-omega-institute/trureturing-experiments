#!/usr/bin/env python3
"""Exact finite checks of a literal forced-ear fan on the original labels."""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path


Chord = tuple[int, int]


def chord(a: int, b: int) -> Chord:
    return (a, b) if a < b else (b, a)


def crosses(e: Chord, f: Chord) -> bool:
    a, b = e
    c, d = f
    return a < c < b < d or c < a < d < b


def independent_sets(n: int, k: int):
    for points in itertools.combinations(range(n), k):
        p = frozenset(points)
        if all((i + 1) % n not in p for i in p):
            yield p


def check_order(n: int) -> dict:
    diagonals = [
        (a, b)
        for a, b in itertools.combinations(range(n), 2)
        if b - a not in (1, n - 1)
    ]
    index = {e: i for i, e in enumerate(diagonals)}
    ear_index = [index[chord((i - 1) % n, (i + 1) % n)] for i in range(n)]
    crossing_masks = [
        sum(1 << j for j, f in enumerate(diagonals) if crosses(e, f))
        for e in diagonals
    ]
    root_stars = [
        sum(1 << j for j, e in enumerate(diagonals) if h in e)
        for h in range(n)
    ]
    complete_mask = (1 << len(diagonals)) - 1
    verified: dict[tuple[frozenset[int], int], tuple[frozenset[int], int]] = {}
    result = {
        "n": n,
        "independent_forced_sets": 0,
        "literal_fans": 0,
        "compatible_root_fans": 0,
        "minimum_compatible_root_crossing_count": n,
        "conditioned_placements": {},
        "positivity_root_witnesses": 0,
        "strict_empty_trace_witnesses": 0,
    }
    for p_size in (2, 3, 4):
        for p in independent_sets(n, p_size):
            result["independent_forced_sets"] += 1
            forced_ears = sum(1 << ear_index[i] for i in p)
            for h in range(n):
                if h in p:
                    continue
                spokes = root_stars[h]
                for j in p:
                    spoke = index.get(chord(h, j))
                    if spoke is not None:
                        spokes &= ~(1 << spoke)
                fan = forced_ears | spokes
                witness = {"n": n, "P": sorted(p), "root": h, "fan_mask": fan}
                assert fan.bit_count() == n - 3, ("diagonal_count", witness)
                assert forced_ears & fan == forced_ears, ("forced", witness)
                for i, crossing in enumerate(crossing_masks):
                    if (fan >> i) & 1:
                        assert not crossing & fan, ("noncrossing", witness, diagonals[i])
                    else:
                        assert crossing & fan, ("maximality", witness, diagonals[i])
                assert fan & ~complete_mask == 0, ("internal_diagonals", witness)
                actual_ears = frozenset(i for i in range(n) if (fan >> ear_index[i]) & 1)
                assert p <= actual_ears, ("actual_forced_ears", witness)
                assert actual_ears <= p | {(h - 1) % n, (h + 1) % n}, (
                    "actual_ear_support", witness, sorted(actual_ears)
                )
                cross_count = (crossing_masks[ear_index[h]] & fan).bit_count()
                if (h - 1) % n not in p and (h + 1) % n not in p:
                    assert cross_count == n - len(p) - 3, ("cross_count", witness)
                    assert cross_count > 1, ("strict_image_exclusion", witness)
                    result["compatible_root_fans"] += 1
                    result["minimum_compatible_root_crossing_count"] = min(
                        result["minimum_compatible_root_crossing_count"], cross_count
                    )
                verified[p, h] = actual_ears, cross_count
                result["literal_fans"] += 1

    for p_size, q_size in ((2, 2), (3, 3), (2, 4), (4, 2)):
        placements = 0
        for k in independent_sets(n, p_size + q_size):
            ek = frozenset(
                h for h in range(n)
                if h not in k and (h - 1) % n not in k and (h + 1) % n not in k
            )
            for points in itertools.combinations(sorted(k), p_size):
                p = frozenset(points)
                q = k - p
                admissible_roots = [
                    h for h in range(n)
                    if h not in p and (h - 1) % n not in q and (h + 1) % n not in q
                ]
                witness = {"n": n, "P": sorted(p), "Q": sorted(q)}
                assert len(p) + 2 * len(q) < n, ("root_budget", witness)
                assert admissible_roots, ("root_supply", witness)
                actual_ears, _ = verified[p, admissible_roots[0]]
                assert p <= actual_ears and actual_ears.isdisjoint(q), (
                    "positive_full_conditioned_fibre", witness
                )
                result["positivity_root_witnesses"] += 1
                for h in ek:
                    actual_ears, cross_count = verified[p, h]
                    potential_trace = ek - {(h - 1) % n, (h + 1) % n}
                    assert h in potential_trace, ("observed_root", witness, h)
                    assert actual_ears.isdisjoint(q), ("strict_forbidden", witness, h)
                    assert actual_ears.isdisjoint(potential_trace), (
                        "uniform_empty_observation", witness, h
                    )
                    assert cross_count > 1, ("outside_flip_image", witness, h)
                    result["strict_empty_trace_witnesses"] += 1
                placements += 1
        result["conditioned_placements"][f"{p_size},{q_size}"] = placements
    return result


def main() -> None:
    assert sys.version_info >= (3, 12), sys.version
    assert __debug__, "Assertions must be enabled"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-n", type=int, default=14)
    parser.add_argument("--max-n", type=int, default=18)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert 14 <= args.min_n <= args.max_n
    rows = [check_order(n) for n in range(args.min_n, args.max_n + 1)]
    result = {
        "schema": "fib-literal-forced-fan-v1",
        "range": [args.min_n, args.max_n],
        "arithmetic": "exact integer labels, unordered internal diagonals and finite bitsets",
        "construction": "ears(P) union all original internal spokes from h to vertices outside P",
        "scope": "all independent P of sizes2,3,4; all source ordered P/Q placements",
        "limits": "finite carrier checks only; no unbounded theorem, flip-existence proof or Lean validation",
        "rows": rows,
        "totals": {
            name: sum(row[name] for row in rows)
            for name in (
                "independent_forced_sets", "literal_fans", "compatible_root_fans",
                "positivity_root_witnesses", "strict_empty_trace_witnesses"
            )
        },
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result["totals"], sort_keys=True))


if __name__ == "__main__":
    main()
