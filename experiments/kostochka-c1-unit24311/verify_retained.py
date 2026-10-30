"""Replay retained C1 evidence with Python 3.13 and the standard library.

Default execution reads fixed inputs and witnesses, replays only the original
small controls, and checks saved table arithmetic. It does not generate graphs
or search for packings on the 30,000 inputs or the large clique sweep.
"""
import argparse
from collections import Counter
import gzip
import hashlib
from itertools import combinations, combinations_with_replacement
import json
from pathlib import Path

from capacity_packing import colorpack
from packing_core import NodeLimit, adjacency, exact_packing, parameters

ROOT = Path(__file__).resolve().parent
EXPECTED_INPUT_SHA = "aebdf7d53ba8177a43a767ec3d4809156fe2ed083c52748b5f28e57696aab875"
EXPECTED_WITNESS_SHA = "3c45d8afd07cc5748a2d694c935b646aff650751151177240a73c521bfec93dd"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def direct_packing(n, edges, sizes, masks):
    """Direct endpoint/cardinality checks, independent of bitset solver code."""
    if len(sizes) != len(masks):
        return False
    used = set()
    for size, mask in zip(sizes, masks):
        if type(mask) is not int or not 0 <= mask < (1 << n):
            return False
        vertices = {v for v in range(n) if mask & (1 << v)}
        if len(vertices) != size or used.intersection(vertices):
            return False
        if any(u in vertices and v in vertices for u, v in edges):
            return False
        used.update(vertices)
    return True


def check_witnesses():
    raw = json.loads((ROOT / "packing-witnesses.json").read_text())
    witnesses = {}
    for seed, index, masks in raw["witnesses"]:
        key = (seed, index)
        require(key not in witnesses, "duplicate witness key")
        witnesses[key] = masks
    input_sha, witness_sha = hashlib.sha256(), hashlib.sha256()
    seen, per_seed = set(), Counter()
    ranges = {key: [None, None] for key in ("ell", "n", "edges", "t")}
    with gzip.open(ROOT / "fixed-inputs.jsonl.gz", "rb") as stream:
        for line in stream:
            row = json.loads(line)
            require(set(row) == {"seed", "index", "ell", "k", "t", "n", "edges"},
                    "input schema")
            require(line == canonical(row), "noncanonical input record")
            seed, index, ell, sizes, t, n, edges = (
                row[key] for key in ("seed", "index", "ell", "k", "t", "n", "edges"))
            require(all(type(x) is int for x in (seed, index, ell, t, n, *sizes)),
                    "integer parameters")
            key = (seed, index)
            require(key not in seen and key in witnesses, "input/witness key mismatch")
            require(seed in (261038430, 1, 2, 3) and 0 <= index < 7500, "fixed input key")
            require(ell in (4, 5, 6) and len(sizes) == ell, "finite ell scope")
            require(sizes[-1] in (2, 3) and all(sizes[-1] <= k <= sizes[-1] + 2
                                              for k in sizes[:-1]), "finite size scope")
            require(0 <= t <= 6 and n == sum(sizes) + t, "n and slack")
            edge_set = set()
            neighbors = [set() for _ in range(n)]
            for edge in edges:
                require(isinstance(edge, list) and len(edge) == 2, "edge shape")
                u, v = edge
                require(type(u) is int and type(v) is int and 0 <= u < v < n,
                        "simple graph endpoint")
                require((u, v) not in edge_set, "duplicate edge")
                edge_set.add((u, v))
                neighbors[u].add(v)
                neighbors[v].add(u)
            require(all(v not in neighbors[v] for v in range(n)), "irreflexivity")
            require(all(v in neighbors[u] for v in range(n) for u in neighbors[v]),
                    "symmetry")
            require(len(edges) < (2 * ell - 1) * (n - sum(sizes) + 1), "strict bound")
            masks = witnesses[key]
            require(direct_packing(n, edges, sizes, masks), f"invalid witness {key}")
            require(direct_packing(n, edges, sizes[:-1], masks[:-1]), "partial witness")
            input_sha.update(line)
            witness_sha.update((json.dumps([seed, index, masks], separators=(",", ":"))
                                + "\n").encode())
            seen.add(key)
            per_seed[str(seed)] += 1
            for name, value in (("ell", ell), ("n", n), ("edges", len(edges)), ("t", t)):
                lo, hi = ranges[name]
                ranges[name] = [value if lo is None else min(lo, value),
                                value if hi is None else max(hi, value)]
    require(len(seen) == len(witnesses) == 30000, "record count")
    require(all(value == 7500 for value in per_seed.values()), "seed count")
    require(input_sha.hexdigest() == EXPECTED_INPUT_SHA, "fixed input stream hash")
    require(witness_sha.hexdigest() == EXPECTED_WITNESS_SHA, "witness stream hash")
    return {"input_count": len(seen), "unique_seed_index_records": len(seen),
            "units": "generated inputs; graph repetitions are allowed",
            "per_seed": dict(sorted(per_seed.items())), "observed_ranges": ranges,
            "valid_full_witnesses": len(seen), "valid_partial_witnesses": len(seen),
            "source_bound_valid": len(seen), "simple_symmetric_irreflexive": len(seen),
            "graph_input_stream_sha256": input_sha.hexdigest(),
            "witness_stream_sha256": witness_sha.hexdigest(),
            "packing_searches_on_saved_inputs": 0}


def naive_exists(n, edges, sizes):
    """Small-control oracle using vertex sets and unordered edge endpoints."""
    choices = {}
    for size in set(sizes):
        choices[size] = []
        for vertices in combinations(range(n), size):
            chosen = frozenset(vertices)
            if not any(u in chosen and v in chosen for u, v in edges):
                choices[size].append(chosen)

    def rec(index, used):
        if index == len(sizes):
            return True
        return any(not used.intersection(chosen) and rec(index + 1, used | chosen)
                   for chosen in choices[sizes[index]])

    return rec(0, frozenset())


def partitions(n, minimum=1):
    if n == 0:
        yield ()
        return
    for size in range(minimum, n + 1):
        for tail in partitions(n - size, size):
            yield (size,) + tail


def clique_edges(parts):
    edges, start = [], 0
    for size in parts:
        edges.extend(combinations(range(start, start + size), 2))
        start += size
    return edges


def clique_prefix(parts, sizes):
    total = 0
    for count, size in enumerate(sorted(sizes, reverse=True), 1):
        total += size
        if total > sum(min(part, count) for part in parts):
            return False
    return True


def clique_all_row_subsets(parts, sizes):
    """Check every row-subset cut; do not rely on the prefix reduction."""
    for count in range(1, len(sizes) + 1):
        for rows in combinations(range(len(sizes)), count):
            if sum(sizes[i] for i in rows) > sum(min(part, count) for part in parts):
                return False
    return True


def check_controls():
    literal_count = color_count = clique_count = 0
    pairs = list(combinations(range(5), 2))
    for bits in range(1 << len(pairs)):
        edges = [edge for i, edge in enumerate(pairs) if bits & (1 << i)]
        a = adjacency(5, edges)
        for sizes in ((2,), (3,), (2, 2), (3, 2)):
            expected = naive_exists(5, edges, sizes)
            literal = exact_packing(a, sizes)
            colored, _ = colorpack(a, sizes)
            require((literal is not None) == expected, "literal/naive disagreement")
            require((colored is not None) == expected, "color/naive disagreement")
            if literal is not None:
                require(direct_packing(5, edges, sorted(sizes, reverse=True), literal),
                        "literal control witness")
            if colored is not None:
                require(direct_packing(5, edges, sizes, colored), "color control witness")
            literal_count += 1
            color_count += 1
    for n in range(2, 11):
        for parts in partitions(n):
            edges = clique_edges(parts)
            a = adjacency(n, edges)
            for ell in range(1, 4):
                for sizes in combinations_with_replacement(range(2, 5), ell):
                    if sum(sizes) > n:
                        continue
                    literal = exact_packing(a, sizes)
                    require(clique_prefix(parts, sizes) == (literal is not None),
                            "clique/literal disagreement")
                    require(clique_all_row_subsets(parts, sizes) == (literal is not None),
                            "all-row-subset/literal disagreement")
                    clique_count += 1
    # Existing fixed-partial obstruction: selecting {2,3} leaves the edge {0,1}.
    edges = [(0, 1)]
    require(direct_packing(4, edges, (2,), (12,)), "fixed partial")
    require(not direct_packing(4, edges, (2,), (3,)), "fixed extension failure")
    for full in (exact_packing(adjacency(4, edges), (2, 2)),
                 colorpack(adjacency(4, edges), (2, 2))[0]):
        require(full is not None and direct_packing(4, edges, (2, 2), full),
                "all-set rechoice witness")
    edges = clique_edges((1, 1, 1, 7))
    a = adjacency(10, edges)
    partial = exact_packing(a, (2, 2, 2))
    require(partial is not None and direct_packing(10, edges, (2, 2, 2), partial),
            "boundary partial")
    require(exact_packing(a, (2, 2, 2, 2)) is None, "boundary literal exclusion")
    require(colorpack(a, (2, 2, 2, 2))[0] is None, "boundary coloring exclusion")
    require(len(edges) == (2 * 4 - 1) * (10 - 8 + 1) == 21, "boundary equality")
    # Existing declared unknown semantics must remain an exception.
    for solver in (exact_packing, lambda a, k, max_nodes: colorpack(a, k, max_nodes)):
        try:
            solver(adjacency(4, []), (2, 2), max_nodes=0)
        except NodeLimit:
            pass
        else:
            raise ValueError("node limit failed to remain unknown")
    require((literal_count, color_count, clique_count) == (4096, 4096, 1751),
            "original control-domain counts")
    return {"literal_vs_endpoint_naive": literal_count,
            "capacity_coloring_vs_endpoint_naive": color_count,
            "clique_prefix_and_all_cuts_vs_literal": clique_count,
            "fixed_partial_fails_but_full_rechoice_succeeds": True,
            "equality_control": {"n": 10, "ell": 4, "k": [2, 2, 2, 2],
                                 "edges": 21, "bound": 21, "partial": True,
                                 "full": False, "strict_source_hypothesis": False},
            "node_limits_are_unknown": True}


def check_table():
    table = json.loads((ROOT / "clique-results.json").read_text())
    rows = table["per_parameter"]
    params = {f"{ell}:{','.join(map(str, sizes))}:{t}": (ell, sizes, t, n)
              for ell, sizes, t, n in parameters()}
    require(rows.keys() == params.keys(), "clique parameter keys")
    # Integer partition counts by dynamic programming; no large graph sweep.
    counts = [1] + [0] * 34
    for part in range(1, 35):
        for n in range(part, 35):
            counts[n] += counts[n - part]
    for n, count in table["partition_counts_by_n"].items():
        require(count == counts[int(n)], "partition count table")
    totals = Counter()
    for key, row in rows.items():
        require(row["all_partitions"] == counts[params[key][3]], "parameter partition count")
        require(0 <= row["strict_edge_eligible"] <= row["all_partitions"], "table bounds")
        totals.update(row)
    saved_keys = set()
    for row in table["minimum_partial_true_full_false"]:
        ell, sizes, t, n, parts = (row[k] for k in ("ell", "k", "t", "n", "parts"))
        key = f"{ell}:{','.join(map(str, sizes))}:{t}"
        require(key in params and key not in saved_keys, "saved clique key")
        require((ell, tuple(sizes), t, n) == params[key], "saved clique parameters")
        require(all(type(c) is int and c > 0 for c in parts) and sum(parts) == n,
                "saved clique components")
        require(row["m"] == sum(c * (c - 1) // 2 for c in parts), "saved edge count")
        require(row["bound"] == (2 * ell - 1) * (t + 1), "saved source bound")
        require(row["gap"] == row["m"] - row["bound"], "saved gap")
        require(clique_all_row_subsets(parts, sizes[:-1]) and
                not clique_all_row_subsets(parts, sizes), "saved cut classification")
        saved_keys.add(key)
    require(saved_keys == params.keys(), "saved clique record count")
    require(totals["all_partitions"] == 778206 and
            totals["strict_edge_eligible"] == 344943, "reported table totals")
    return {"parameter_rows": len(rows), "table_all_partition_pairs_sum": 778206,
            "table_strict_eligible_sum": 344943,
            "partition_counts_dp_checked_n": [8, 34],
            "saved_nonpacking_clique_records_cut_checked": len(saved_keys),
            "large_clique_sweep_reexecuted": False,
            "unverified": ["actual sweep traversal", "all strict-eligible classifications",
                           "claimed minima over partitions", "arbitrary-graph coverage"]}


def check_upstream_bytes():
    expected = json.loads((ROOT / "upstream-artifact-sha256.json").read_text())
    checked = {}
    for name in ("packing_core.py", "packing-witnesses.json", "clique-results.json"):
        data = (ROOT / name).read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        require(len(data) == expected[name]["bytes"] and digest == expected[name]["sha256"],
                f"original artifact bytes changed: {name}")
        checked[name] = digest
    return checked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write a deterministic validation report")
    parser.add_argument("--check-report", type=Path, help="compare against a saved report")
    args = parser.parse_args()
    report = {"whole_answer_status": "unresolved", "candidate": None,
              "original_artifact_file_sha256": check_upstream_bytes(),
              "witnesses": check_witnesses(), "controls": check_controls(),
              "saved_clique_table": check_table(),
              "limits": ["finite evidence only", "no new search or solver route",
                         "no source refetch or new literature/library search",
                         "ell boundary remains unresolved for positive whole-source closure",
                         "no Lean/kernel closure, admission, completion, novelty or KPI"]}
    if args.check_report:
        require(json.loads(args.check_report.read_text()) == report, "saved report mismatch")
    if args.output:
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
