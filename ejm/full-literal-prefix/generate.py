"""Propose the exact original 201-context suffix with a full literal cache.

The input layout is the existing actual-incremental-cache extraction consumed
by check-original-incremental-cache.cpp. Generation supplies data, never a Lean
receipt: the emitted full-field equality and unchanged suffix require kernel
checking. This program does not invoke a compiler or modify existing modules.
"""

import argparse
import hashlib
from pathlib import Path
import re


def matrix(tokens, height, width):
    return [[int(next(tokens)) for _ in range(width)] for _ in range(height)]


def read_inputs(directory):
    with (directory / "actual-incremental-cache.txt").open() as stream:
        tokens = (word for line in stream for word in line.split())
        assert next(tokens) == "EJM_ACTUAL_INCREMENTAL_CACHE_V1"
        outcomes = matrix(tokens, 24, 4)
        outcome_mul = matrix(tokens, 24, 24)
        inverse = [int(next(tokens)) for _ in range(12)]
        blocks = []
        for expected_n, expected_count in [(10, 41294), (9, 8689)]:
            n, count = int(next(tokens)), int(next(tokens))
            assert (n, count) == (expected_n, expected_count)
            permutations = matrix(tokens, 12, n)
            coordinate_mul = matrix(tokens, 12, 12)
            incidence = matrix(tokens, 20, 16)
            rows = []
            for _ in range(count):
                key, lower = int(next(tokens)), int(next(tokens))
                rows.append((key, lower, matrix(tokens, n, 4)))
            blocks.append(dict(n=n, count=count, permutations=permutations,
                               coordinate_mul=coordinate_mul,
                               incidence=incidence, rows=rows))
        assert next(tokens, None) is None
    return outcomes, outcome_mul, inverse, blocks


def read_mirrors(source, label, block):
    text = (source / f"EJM{label}CacheData.lean").read_text()
    matches = re.findall(
        r"\(\.node (\d+) ⟨(\d+), \((-?\d+)\), (0x[0-9a-f]+), (0x[0-9a-f]+)⟩",
        text)
    mirrors = {int(j): (int(key), int(lower), int(slots, 16), int(edges, 16))
               for j, key, lower, slots, edges in matches}
    assert len(matches) == block["count"]
    assert sorted(mirrors) == list(range(block["count"]))
    for index, (key, lower, words) in enumerate(block["rows"]):
        mkey, mlower, slots, edges = mirrors[index]
        assert (mkey, mlower) == (key, lower)
        for coordinate, values in enumerate(words):
            slot = slots >> (4 * coordinate) & 15
            for value, word in enumerate(values):
                assert word == edges >> (25 * (4 * slot + value)) & ((1 << 25) - 1)
    return mirrors


def states_after_zero_prefix(outcomes, outcome_mul, inverse, blocks, mirrors):
    states = [[(block["count"] - 1, 0, 0) for _ in range(16)]
              for block in blocks]
    for coordinate in range(15):
        for b, block in enumerate(blocks):
            for k, state in enumerate(states[b]):
                local_coordinate = block["incidence"][coordinate][k]
                if local_coordinate < 0:
                    continue
                index, c, o = state
                transported = block["permutations"][inverse[c]][local_coordinate]
                key, _, words = block["rows"][index]
                assert key // 5 ** transported % 5 == 4
                word = words[transported][outcomes[o][0]]
                next_index = word % 65536
                action = word // 65536 % 512
                assert next_index < block["count"] and action < 288
                states[b][k] = (next_index,
                                block["coordinate_mul"][c][action // 24],
                                outcome_mul[action % 24][o])
    # Check each complete transported local partial, in addition to row IDs.
    for b, block in enumerate(blocks):
        for k, (index, c, o) in enumerate(states[b]):
            local = [None] * block["n"]
            for global_coordinate in range(20):
                i = block["incidence"][global_coordinate][k]
                if i >= 0:
                    assert local[i] is None
                    local[i] = 0 if global_coordinate < 15 else 4
            assert all(value is not None for value in local)
            transported = [4 if local[i] == 4 else outcomes[o][local[i]]
                           for i in block["permutations"][c]]
            key = mirrors[b][index][0]
            assert transported == [key // 5 ** i % 5 for i in range(block["n"])]
    return states


def source_text(states, mirrors):
    vectors = []
    for b in range(2):
        values = []
        for index, c, o in states[b]:
            key, lower, slots, edges = mirrors[b][index]
            values.append(f"    ⟨({index}, ({c}, {o})), "
                          f"⟨{key}, ({lower}), 0x{slots:x}, 0x{edges:x}⟩⟩")
        vectors.append("![\n" + ",\n".join(values) + "]")
    literal = "⟨" + ",\n  ".join(vectors) + "⟩"
    return '''/- Prepared cost comparison, not yet kernel checked. The complete
literal cache must first equal the actual fifteen-step prefix. Both the cache
equality and the unchanged original 201-context suffix require Lean receipts.
This sample asserts no complete global lower bound or locality exclusion. -/
import D5.S3.Quantum.Entanglement.EJM.EJMOriginalCachedGlobalReplay
import D5.S3.Quantum.Entanglement.EJM.EJMCacheMaterialization

set_option autoImplicit false
set_option relaxedAutoImplicit false
set_option maxHeartbeats 100000000
set_option maxRecDepth 65536
set_option linter.style.longLine false

namespace D5.S3.Quantum.Entanglement.EJM.EJMOriginalGlobalSuffixMaterializedPilot
open D5.S0.Certificates.FiniteSymmetryBoundReplay (Partial)
open D5.S0.Certificates.FiniteCachedSymmetryDAGReplay (check)
open D5.S3.Quantum.Entanglement.EJM.EJMCacheState
open D5.S3.Quantum.Entanglement.EJM.EJMCacheSoundness (initial)
open D5.S3.Quantum.Entanglement.EJM.EJMCacheMaterialization (checkSame checked_same)
open D5.S3.Quantum.Entanglement.EJM.EJMPartialBoundData (scale)

private def prefix : Fin 15 → Fin 4 := ![0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

def computedCache : Cache :=
  (List.finRange 15).foldl
    (fun cache i => step cache ⟨i.val, by have h := i.isLt; omega⟩ (prefix i))
    initial

def partialEvent : Partial (Fin 20) (Fin 4) := fun g =>
  if h : g.val < 15 then some (prefix ⟨g.val, h⟩) else none

/-- Every state and row payload is literal; no mirror lookup is left here. -/
def literalCache : Cache := ''' + literal + '''

set_option trace.profiler true in
private theorem full_fields : checkSame literalCache computedCache = true := by
  decide +kernel

/-- Full equality includes all thirty-two states and every future row field. -/
theorem literal_cache_exact : literalCache = computedCache :=
  checked_same _ _ full_fields

set_option trace.profiler true in
theorem literal_suffix :
    check EJMOriginalGlobalDAGData.nodes step lower EJMGlobalScoreSymmetry.dominates
      scale 6 ⟨55, by decide⟩ literalCache partialEvent = true := by
  decide +kernel

/-- The result returns to the actual original initial cache and all fifteen updates. -/
theorem actual_original_context :
    check EJMOriginalGlobalDAGData.nodes step lower EJMGlobalScoreSymmetry.dominates
      scale 6 ⟨55, by decide⟩ computedCache partialEvent = true := by
  simpa only [literal_cache_exact] using literal_suffix

#print axioms literal_cache_exact
#print axioms literal_suffix
#print axioms actual_original_context
end D5.S3.Quantum.Entanglement.EJM.EJMOriginalGlobalSuffixMaterializedPilot
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    outcomes, outcome_mul, inverse, blocks = read_inputs(args.inputs)
    mirrors = [read_mirrors(args.source_dir, label, block)
               for label, block in zip(["Triangle", "Singleton"], blocks)]
    states = states_after_zero_prefix(outcomes, outcome_mul, inverse, blocks, mirrors)
    text = source_text(states, mirrors)
    assert "mirror " not in text.split("def literalCache : Cache :=", 1)[1].split(
        "set_option trace.profiler", 1)[0]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text)
    print(args.output, "sha256=" + hashlib.sha256(text.encode()).hexdigest(),
          "bytes=" + str(len(text.encode())), "kernel_accepted=false")


if __name__ == "__main__":
    main()
