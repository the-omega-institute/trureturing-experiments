"""Bounded falsifiers for packed_constructor.py; no evaluator execution."""
from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path

from packed_constructor import SHAPES, PackedState

assert sys.version_info >= (3, 12) and __debug__


def run() -> dict:
    started = time.perf_counter()
    state = PackedState(3)
    shape_handles = []
    for schema, shape in enumerate(SHAPES):
        handle = state.allocate(0, shape, schema)
        shape_handles.append(handle)
        for bit in range(shape[1]):
            state.write_bit(handle, bit, bit & 1)
    alias = state.allocate(1, (4, 8), 11)
    cycle = state.allocate(2, (2, 0), 12)
    state.write_ref(alias, 0, shape_handles[0])
    state.write_ref(alias, 1, cycle)
    state.write_ref(alias, 2, cycle)
    state.write_ref(cycle, 0, alias)
    state.header.pc = alias
    state.header.old_pc = shape_handles[0]
    state.header.pending = cycle
    state.header.cursor = alias
    free = state.allocate(0, (1, 0), 13)
    state.release_to_free(free)
    before = state.snapshot()
    assert state.read_ref(alias, 1) == state.read_ref(alias, 2) == cycle
    assert state.read_ref(cycle, 0) == alias
    assert state.read_bit(alias, 1) == 0
    state.write_bit(alias, 1, 1)
    assert state.read_bit(alias, 1) == 1
    source = [state.allocate(1, (4, 8), 20) for _ in range(3)]
    target = [state.allocate(1, (4, 8), 21) for _ in range(3)]
    for i, handle in enumerate(source):
        state.write_ref(handle, 0, source[(i + 1) % 3])
        state.write_bit(handle, 0, i & 1)
    state.relocate(source, target)
    for old, new in zip(source, target):
        assert state.read_ref(new, 0) == target[(source.index(state.read_ref(old, 0)))]
        assert state.read_bit(new, 0) == state.read_bit(old, 0)
    mutated = state.snapshot()
    state.write_ref(alias, 2, shape_handles[1])
    assert state.snapshot() != mutated
    state.write_ref(alias, 2, cycle)
    state.width.update({0: 4, 1: 8, 2: 12})
    census = state.modeled_census()
    assert all(row["objects"] == row["grant_bits"] for row in census.values())
    assert before != state.snapshot()
    return {
        "schema": 1,
        "shapes": len(SHAPES),
        "slots": sum(len(slots) for slots in state.owners.values()),
        "alias_cycle_pending_free": True,
        "installation_forward_aliases": True,
        "alias_mutation_detected": True,
        "owner_widths": state.width,
        "census": census,
        "logical_packed_bits": state.packed_bits(),
        "host_rss_ru_maxrss": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "wall_seconds": time.perf_counter() - started,
        "events": len(state.events),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = json.dumps(run(), sort_keys=True, indent=2) + "\n"
    if args.output:
        tmp = args.output.with_suffix(args.output.suffix + ".tmp")
        tmp.write_text(result)
        tmp.replace(args.output)
    print(result, end="")


if __name__ == "__main__":
    main()
