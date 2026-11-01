"""Finite checks of the actual diagonal-set map in B-37.2, equation (37.8).

Enumerate triangulations independently of the slide map and compare its image
with the full independently filtered target. These finite checks do not prove
the unbounded theorem or certify Lean source.

Python >= 3.12; standard library only.
"""

import argparse
from functools import cache
from itertools import combinations
import json
from math import comb
from pathlib import Path
import sys

assert sys.version_info >= (3, 12)
Edge = tuple[int, int]
Triangulation = frozenset[Edge]


def edge(a: int, b: int) -> Edge:
    return min(a, b), max(a, b)


def diagonal(e: Edge, n: int) -> bool:
    a, b = e
    return a < b and b - a not in (1, n - 1)


def crosses(e: Edge, f: Edge) -> bool:
    a, b = e
    c, d = f
    return a < c < b < d or c < a < d < b


@cache
def interval(a: int, b: int) -> frozenset[Triangulation]:
    if b <= a + 2:
        return frozenset({frozenset()})
    out = set()
    for k in range(a + 1, b):
        seams = set()
        if k > a + 1:
            seams.add(edge(a, k))
        if b > k + 1:
            seams.add(edge(k, b))
        for left in interval(a, k):
            for right in interval(k, b):
                out.add(left | right | seams)
    return frozenset(out)


def ear(k: int, n: int) -> Edge:
    return edge((k - 1) % n, (k + 1) % n)


def swap_edge(e: Edge, x: int, y: int) -> Edge:
    def swap(v: int) -> int:
        return y if v == x else x if v == y else v
    return edge(swap(e[0]), swap(e[1]))


def slide(t: Triangulation, n: int, start: int, inverse: bool = False) -> Triangulation:
    x, y = (start + 1) % n, (start + 2) % n
    s, old_t = ear(x, n), ear(y, n)
    erased, inserted = (s, old_t) if inverse else (old_t, s)
    return frozenset(swap_edge(e, x, y) for e in t - {erased}) | {inserted}


def check_enumeration(n: int) -> dict:
    all_t = interval(0, n - 1)
    assert len(all_t) == comb(2 * (n - 2), n - 2) // (n - 1)
    diagonals = {e for e in combinations(range(n), 2) if diagonal(e, n)}
    for t in all_t:
        assert len(t) == n - 3 and t <= diagonals
        assert all(not crosses(e, f) for e, f in combinations(t, 2))
        assert all(any(crosses(e, f) for f in t) for e in diagonals - t)
    applications = 0
    for start in range(n):
        s, old_t = ear((start + 1) % n, n), ear((start + 2) % n, n)
        source = {t for t in all_t if old_t in t}
        target = {t for t in all_t if s in t}
        image = set()
        remote = set(range(n)) - {(start + k) % n for k in range(4)}
        for t in source:
            u = slide(t, n, start)
            assert u in target, (n, start, t, u)
            assert slide(u, n, start, inverse=True) == t
            assert all((ear(k, n) in t) == (ear(k, n) in u) for k in remote)
            assert ear((start + 1) % n, n) not in t
            assert ear((start + 3) % n, n) not in t
            assert ear(start, n) not in u
            assert ear((start + 2) % n, n) not in u
            image.add(u)
        assert image == target, (n, start, "full-target mismatch")
        assert all(slide(slide(u, n, start, inverse=True), n, start) == u for u in target)
        applications += len(source)
    return {"n": n, "triangulations": len(all_t), "block_placements": n,
            "source_applications": applications}


def check_local_transport(n: int) -> dict:
    chords = [e for e in combinations(range(n), 2) if 2 not in e]
    mapped = {e: swap_edge(e, 1, 2) for e in chords}
    assert all(crosses(e, f) == crosses(mapped[e], mapped[f])
               for e, f in combinations(chords, 2))
    assert all(not diagonal(e, n) or e == (1, 3) or diagonal(mapped[e], n)
               for e in chords)
    for e in combinations(range(n), 2):
        if diagonal(e, n) and 2 in e:
            assert crosses(e, (1, 3))
    return {"n": n, "chords_avoiding_deleted_vertex": len(chords),
            "crossing_pairs": comb(len(chords), 2)}


def forced_completions(n: int, forced: set[int]):
    """Root-triangle enumeration on the polygon after independent ear deletion."""
    assert all((a - b) % n not in (1, n - 1) for a, b in combinations(forced, 2))
    vertices = [v for v in range(n) if v not in forced]
    inserted = frozenset(ear(v, n) for v in forced)
    for reduced in interval(0, len(vertices) - 1):
        yield frozenset(edge(vertices[a], vertices[b]) for a, b in reduced) | inserted


def check_conditioned(n: int, p: int, q: int, keep_u: bool) -> dict:
    """Check every exact observed trace for one literal compatible P/Q placement."""
    start = 5 if keep_u else 4
    k = set(range(start, start + 2 * (p + q), 2))
    assert max(k) < n and p + q <= n // 2
    forced = set(sorted(k)[:p])
    forbidden = k - forced
    k_edges = {ear(v, n) for v in k}
    common = {v for v in range(n) if v not in k and
              all(not crosses(ear(v, n), e) for e in k_edges)}
    observed = common - {0}
    assert {1, 2} <= observed
    assert (3 in observed) == keep_u
    remote_observed = observed - {1, 2, 3}
    forbidden_edges = {ear(v, n) for v in forbidden}
    source_groups: dict[tuple[int, ...], set[Triangulation]] = {}
    target_groups: dict[tuple[int, ...], set[Triangulation]] = {}
    counts = []
    for centre, groups in [(2, source_groups), (1, target_groups)]:
        count = 0
        for t in forced_completions(n, forced | {centre}):
            if t & forbidden_edges:
                continue
            assert len(t) == n - 3
            assert all(diagonal(e, n) for e in t)
            assert all(not crosses(e, f) for e, f in combinations(t, 2))
            trace = {v for v in range(n) if ear(v, n) in t}
            r = tuple(sorted(trace & remote_observed))
            if centre == 2:
                assert trace & observed == set(r) | {2}
                image = slide(t, n, 0)
                assert slide(image, n, 0, inverse=True) == t
                groups.setdefault(r, set()).add(image)
            else:
                assert trace & (observed - {2, 3}) == set(r) | {1}
                assert slide(slide(t, n, 0, inverse=True), n, 0) == t
                groups.setdefault(r, set()).add(t)
            count += 1
        counts.append(count)
    assert source_groups == target_groups, (n, p, q, keep_u)
    assert counts[0] == counts[1]
    return {"n": n, "p": p, "q": q, "u_observed": keep_u,
            "P": sorted(forced), "Q": sorted(forbidden), "J": sorted(observed),
            "nonempty_trace_fibres": len(source_groups),
            "source_completions": counts[0], "target_completions": counts[1]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-enumeration-n", type=int, default=11)
    parser.add_argument("--max-local-n", type=int, default=30)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not 5 <= args.max_enumeration_n <= 12 or args.max_local_n < 5:
        parser.error("enumeration range must end at 5..12; local range must reach 5")
    enumeration = [check_enumeration(n) for n in range(5, args.max_enumeration_n + 1)]
    local = [check_local_transport(n) for n in range(5, args.max_local_n + 1)]
    pentagon = [frozenset({(1, 3), (0, 3)}), frozenset({(1, 3), (1, 4)})]
    expected = [frozenset({(0, 2), (0, 3)}), frozenset({(0, 2), (2, 4)})]
    assert [slide(t, 5, 0) for t in pentagon] == expected
    conditioned = [check_conditioned(14, 2, 2, True),
                   check_conditioned(14, 2, 2, False)]
    for p, q in [(3, 3), (2, 4), (4, 2)]:
        conditioned += [check_conditioned(15, p, q, False),
                        check_conditioned(16, p, q, True)]
    result = {"enumeration": enumeration, "local_transport": local,
              "conditioned_fibres": conditioned,
              "all_assertions_passed": True,
              "limits": ["Finite experiment; no unbounded or Lean proof.",
                         "Conditioned checks use eight literal P/Q/J placements; the general source-condition discharge remains a proof obligation."]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"all_assertions_passed": True,
                      "source_applications": sum(x["source_applications"] for x in enumeration),
                      "crossing_pairs": sum(x["crossing_pairs"] for x in local),
                      "conditioned_placements": len(conditioned),
                      "conditioned_source_completions": sum(x["source_completions"] for x in conditioned),
                      "enumeration_range": [5, args.max_enumeration_n],
                      "local_range": [5, args.max_local_n]}))


if __name__ == "__main__":
    main()
