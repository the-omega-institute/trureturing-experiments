"""Bounded source-structure diagnostic, never Lean acceptance or timing.

Read the unchanged original 55/15/0 suffix (201 actual contexts), not the
87M-context complete expansion. Counts distinguish source function chains
from an eager vector schedule. Kernel reduction sharing is deliberately not
assumed, so neither count is a prediction of kernel operations or speedup.
The DAG is memory mapped and only bounded records are decoded.
"""
import argparse
from collections import Counter
import hashlib
import json
import mmap
from pathlib import Path
import re
import struct


DAG_SHA = "bf188b5ef8d359fc12872297ba3c3db9b633debdec7db6eddf924265c880bf59"
DEFAULT_DATA = Path("/Users/macstudio/omega/.formalization-goal-01a12631/checkpoints/"
                    "ejm-mdd-global/existing-cover-original-suffix-alignment-20261014/"
                    "original-source-identity-stream")
DEFAULT_SOURCE = Path("/Users/macstudio/omega/trureturing-01a12631-certificate/"
                      "D5/S3/Quantum/Entanglement/EJM")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=DEFAULT_DATA)
    ap.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    cache_source = (args.source / "EJMLocalCacheData.lean").read_text()
    incidence = []
    for name, following in [("triangleIncidence", "def singletonIncidence"),
                            ("singletonIncidence", "end D5.")]:
        body = cache_source.split("def " + name, 1)[1].split(following, 1)[0]
        entries = re.findall(r"some (\d+)|\b(none)\b", body)
        values = [int(some) if some else -1 for some, _ in entries]
        assert len(values) == 20 * 16
        incidence.append([values[16*g:16*g+16] for g in range(20)])
    affected = [sum(i >= 0 for block in incidence for i in block[g]) for g in range(20)]

    file = (args.data / "dag-87458265.bin").open("rb")
    dag = mmap.mmap(file.fileno(), 0, access=mmap.ACCESS_READ)
    assert struct.unpack_from("<II", dag) == (0x454A4D32, 1558777)
    assert len(dag) == 8 + 24 * 1558777
    # Full-byte identity is fast streaming, not full-context traversal.
    assert hashlib.sha256(dag).hexdigest() == DAG_SHA
    kinds, depths, indices = Counter(), Counter(), Counter()
    bounds = []
    splits = []

    def walk(index, depth, code):
        assert 15 <= depth <= 20 and 0 <= code < 4 ** depth
        kind, coordinate, *children = struct.unpack_from("<6I", dag, 8 + 24 * index)
        assert kind in (0, 1, 2)
        kinds[kind] += 1
        depths[(kind, depth)] += 1
        indices[index] += 1
        if kind == 0:
            bounds.append((depth, code, index))
        elif kind == 2:
            assert coordinate == depth and all(c < index for c in children)
            splits.append((depth, code, index))
            for value, child in enumerate(children):
                walk(child, depth + 1, 4 * code + value)

    walk(55, 15, 0)
    assert sum(kinds.values()) == 201
    assert 1 + 4 * kinds[2] == 201
    full_chain_incidence = 32 * sum(depth for depth, _, _ in bounds)
    suffix_chain_incidence = 32 * sum(depth - 15 for depth, _, _ in bounds)
    full_chain_transitions = sum(sum(affected[:depth]) for depth, _, _ in bounds)
    suffix_chain_transitions = sum(sum(affected[15:depth]) for depth, _, _ in bounds)
    eager_suffix_incidence = 32 * 4 * len(splits)
    eager_suffix_transitions = 4 * sum(affected[depth] for depth, _, _ in splits)
    result = dict(
        status="bounded_source_structure_only_not_kernel_or_runtime_profile",
        context=dict(index=55, depth=15, code=0, fuel=6, contexts=201),
        dag_sha256=DAG_SHA,
        source_sha256={name: hashlib.sha256((args.source / name).read_bytes()).hexdigest()
                       for name in ["EJMCacheState.lean", "EJMLocalCacheData.lean",
                                    "EJMOriginalGlobalDAGData.lean"]},
        kinds={"bound": kinds[0], "symmetry": kinds[1], "split": kinds[2]},
        distinct_shapes=len(indices),
        actual_contexts_per_shape={str(j): n for j, n in sorted(indices.items())},
        bound_depth_histogram={str(d): n for (k, d), n in sorted(depths.items()) if k == 0},
        split_depth_histogram={str(d): n for (k, d), n in sorted(depths.items()) if k == 2},
        affected_blocks_per_coordinate=affected,
        bound_state_reads=32 * len(bounds),
        function_chain_expansion_without_sharing=dict(
            full_prefix_and_suffix_incidence_tests=full_chain_incidence,
            suffix_only_incidence_tests=suffix_chain_incidence,
            full_prefix_and_suffix_local_transitions=full_chain_transitions,
            suffix_only_local_transitions=suffix_chain_transitions),
        eager_vector_schedule_all_children=dict(
            prefix_incidence_tests=32*15,
            prefix_local_transitions=sum(affected[:15]),
            suffix_incidence_tests=eager_suffix_incidence,
            suffix_local_transitions=eager_suffix_transitions,
            bound_state_reads=32*len(bounds)),
        limits=[
            "No kernel memoization, whnf sharing, vector eagerness or projection cost was measured.",
            "A symmetry leaf does not demand its cache; the eager-all-children count can overstate vector work.",
            "Function-chain counts count one projected-state reconstruction per bound field, not duplication inside triangleStep/singletonStep.",
            "These structural counts justify a controlled representation experiment, not a speedup claim.",
            "All original node occurrences and four value branches are traversed; no acceptance is shared by shape."])
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    dag.close()
    file.close()
    print(json.dumps({k: result[k] for k in ["context", "kinds", "distinct_shapes",
          "bound_state_reads", "function_chain_expansion_without_sharing",
          "eager_vector_schedule_all_children"]}))


if __name__ == "__main__":
    main()
