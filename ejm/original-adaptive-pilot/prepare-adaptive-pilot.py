"""Prepare adaptive actual-context receipts for the unchanged original EJM DAG.

The default context is the complete root 1558776 at depth/code 0/0. Each
suffix retains the entire original prefixCache, prefixPartial, fuel and index.
Only literal node equalities share a (depth, index) key; Acceptance never does.
All four children are emitted with codes 4*code+a, including repeated shapes.
Generated modules remain prepared/uncompiled until canonical kernel checking.
No compiler is invoked. --output-root is required; regression uses 55/15/0.
"""

import argparse
from array import array
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
import struct


DAG_SHA256 = "bf188b5ef8d359fc12872297ba3c3db9b633debdec7db6eddf924265c880bf59"
SIZE = 1558777
GLOBAL_ROOT = SIZE - 1
GLOBAL_CONTEXTS = 87458265
NS = "D5.S3.Quantum.Entanglement.EJM"
MODULE_PATH = Path(*NS.split("."))
DEFAULT_DAG = Path(
    "/Users/macstudio/omega/.formalization-goal-01a12631/checkpoints/ejm-mdd-global/"
    "existing-cover-original-suffix-alignment-20261014/"
    "original-source-identity-stream/dag-87458265.bin"
)


class Dag:
    def __init__(self, path):
        self.raw = path.read_bytes()
        assert hashlib.sha256(self.raw).hexdigest() == DAG_SHA256
        assert struct.unpack_from("<II", self.raw) == (0x454A4D32, SIZE)
        assert len(self.raw) == 8 + 24 * SIZE
        # Count expanded contexts, not distinct shapes. Repeated children
        # contribute once per original value branch. No recursive recounting.
        self.counts = array("Q", [0]) * SIZE
        for index in range(SIZE):
            kind, coordinate, *children = self.row(index)
            assert kind in (0, 1, 2)
            if kind == 2:
                assert 0 <= coordinate < 20 and all(c < index for c in children)
                self.counts[index] = 1 + sum(self.counts[c] for c in children)
            else:
                self.counts[index] = 1
        assert self.counts[55] == 201
        assert self.counts[GLOBAL_ROOT] == GLOBAL_CONTEXTS

    def row(self, index):
        assert 0 <= index < SIZE
        return struct.unpack_from("<6I", self.raw, 8 + 24 * index)


@dataclass(frozen=True)
class Ref:
    name: str
    part: int


def fin(index):
    return f"⟨{index}, by decide⟩"


class Parts:
    """Keep only one bounded part in memory; imports follow actual proof uses."""

    def __init__(self, output_root, name, limit):
        self.directory = output_root / MODULE_PATH
        self.directory.mkdir(parents=True, exist_ok=True)
        self.name, self.limit = name, limit
        self.part = 0
        self.bodies, self.dependencies = [], set()
        self.declarations = 0
        self.targets = (output_root / "targets.txt").open("w")
        self.shards = (output_root / "shards.jsonl").open("w")

    def module(self, part):
        return f"{NS}.{self.name}Part{part:06d}"

    def header(self, imports):
        return [
            "/- Generated actual-context certificate source; prepared/uncompiled.",
            "No native count, minimum or producer assertion is a Lean premise. -/",
            f"import {NS}.EJMOriginalReplayAssembly",
            *(f"import {module}" for module in imports),
            "set_option autoImplicit false",
            "set_option relaxedAutoImplicit false",
            "set_option maxHeartbeats 100000000",
            "set_option maxRecDepth 65536",
            f"namespace {NS}.{self.name}",
            f"open {NS}.EJMOriginalReplayAssembly",
            "open D5.S0.Certificates.FiniteCachedSymmetryDAGReplay (checked_split_node)",
            "",
        ]

    def add(self, name, lines, requires=()):
        if len(self.bodies) == self.limit:
            self.flush()
        self.dependencies.update(ref.part for ref in requires if ref.part != self.part)
        self.bodies.append("\n".join(lines))
        self.declarations += 1
        return Ref(name, self.part)

    def flush(self):
        if not self.bodies:
            return
        assert all(part < self.part for part in self.dependencies)
        imports = [self.module(part) for part in sorted(self.dependencies)]
        text = "\n".join(self.header(imports) + self.bodies + [
            f"end {NS}.{self.name}", ""
        ])
        path = self.directory / f"{self.name}Part{self.part:06d}.lean"
        path.write_text(text)
        self.targets.write(self.module(self.part) + "\n")
        self.shards.write(json.dumps(dict(
            module=self.module(self.part), declarations=len(self.bodies),
            imports=imports, sha256=hashlib.sha256(text.encode()).hexdigest(),
        )) + "\n")
        self.part += 1
        self.bodies, self.dependencies = [], set()

    def finish(self, root, depth, code, index):
        self.flush()
        module = f"{NS}.{self.name}Root"
        lines = self.header([self.module(root.part)])
        if (depth, code, index) == (0, 0, GLOBAL_ROOT):
            lines += [
                "theorem original_globalAccepted :",
                "    EJMOriginalCachedGlobalReplay.globalAccepted := by",
                f"  change Accepted 0 0 {fin(GLOBAL_ROOT)}",
                f"  exact {root.name}", "",
                "#print axioms original_globalAccepted",
            ]
        else:
            lines += [
                f"theorem accepted_actual_context : Accepted {depth} {code} {fin(index)} :=",
                f"  {root.name}", "", "#print axioms accepted_actual_context",
            ]
        lines += [f"end {NS}.{self.name}", ""]
        path = self.directory / f"{self.name}Root.lean"
        path.write_text("\n".join(lines))
        self.targets.write(module + "\n")
        self.targets.close()
        self.shards.close()
        return module, hashlib.sha256(path.read_bytes()).hexdigest()


def generate(dag, args):
    assert 0 <= args.root < SIZE and 0 <= args.depth <= 20
    assert 0 <= args.code < 4 ** args.depth
    assert args.cap >= 1 and args.declarations_per_part >= 1
    assert re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", args.name)
    # Refuse to mingle a new certificate with existing generated proof sources.
    args.output_root.mkdir(parents=True, exist_ok=False)
    parts = Parts(args.output_root, args.name, args.declarations_per_part)
    node_proofs = {}
    suffixes = splits = 0
    represented_contexts = 0

    with (args.output_root / "contexts.jsonl").open("w") as contexts:
        def emit(index, depth, code):
            nonlocal suffixes, splits, represented_contexts
            label = f"accepted_{depth}_{code}_{index}"
            count = dag.counts[index]
            kind, coordinate, *children = dag.row(index)
            if count <= args.cap:
                suffixes += 1
                represented_contexts += count
                # Spell the complete unchanged cache/partial check at each cut.
                ref = parts.add(label, [
                    f"theorem {label} : Accepted {depth} {code} {fin(index)} := by",
                    "  change D5.S0.Certificates.FiniteCachedSymmetryDAGReplay.check",
                    "    EJMOriginalGlobalDAGData.nodes EJMCacheState.step EJMCacheState.lower",
                    "    EJMGlobalScoreSymmetry.dominates EJMPartialBoundData.scale",
                    f"    {21-depth} {fin(index)} (prefixCache {depth} {code})",
                    f"    (prefixPartial {depth} {code}) = true",
                    "  decide +kernel", "",
                ])
                role = "suffix"
            else:
                assert kind == 2 and coordinate == depth and depth < 20
                splits += 1
                represented_contexts += 1
                # Deliberately never memoize emit/Accepted by shape or index.
                branches = [emit(child, depth + 1, 4 * code + value)
                            for value, child in enumerate(children)]
                key = (depth, index)
                vector = "![" + ",".join(fin(child) for child in children) + "]"
                if key not in node_proofs:
                    node_label = f"node_{depth}_{index}"
                    node_proofs[key] = parts.add(node_label, [
                        f"theorem {node_label} :",
                        f"    EJMOriginalGlobalDAGData.nodes {21-depth} {fin(index)} =",
                        f"      .split (coordinate {depth}) ({vector}) := by",
                        "  apply checked_split_node", "  decide +kernel", "",
                    ])
                node = node_proofs[key]
                ref = parts.add(label, [
                    f"theorem {label} : Accepted {depth} {code} {fin(index)} := by",
                    f"  apply assemble {depth} {code} {fin(index)} (by decide) ({vector}) {node.name}",
                    "  intro a", "  fin_cases a",
                    *("  · exact " + child.name for child in branches), "",
                ], requires=[node, *branches])
                role = "split"
            contexts.write(json.dumps(dict(
                index=index, depth=depth, code=code, fuel=21-depth, contexts=count,
                role=role, theorem=label, module=parts.module(ref.part),
                prefix_cache=f"prefixCache {depth} {code}",
                prefix_partial=f"prefixPartial {depth} {code}",
                children=[dict(index=child, depth=depth+1, code=4*code+value)
                          for value, child in enumerate(children)] if role == "split" else [],
            )) + "\n")
            return ref

        root = emit(args.root, args.depth, args.code)
    assert suffixes == 1 + 3 * splits
    assert represented_contexts == dag.counts[args.root]
    target, sha256 = parts.finish(root, args.depth, args.code, args.root)
    record = dict(
        status="prepared/uncompiled actual-context adaptive certificate source",
        kernel_accepted=False, full_global_accepted=False,
        dag_sha256=DAG_SHA256, shapes=SIZE, global_root=GLOBAL_ROOT,
        global_expanded_contexts=dag.counts[GLOBAL_ROOT],
        root=args.root, depth=args.depth, code=args.code,
        full_contexts=dag.counts[args.root], suffix_cap=args.cap,
        suffixes=suffixes, splits=splits, node_equalities=len(node_proofs),
        declarations_per_part=args.declarations_per_part,
        declarations=parts.declarations, parts=parts.part,
        root_target=target, root_source_sha256=sha256,
        emits_original_globalAccepted=(args.depth, args.code, args.root) == (0, 0, GLOBAL_ROOT),
    )
    (args.output_root / "manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dag", type=Path, default=DEFAULT_DAG)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--root", type=int, default=GLOBAL_ROOT)
    parser.add_argument("--depth", type=int, default=0)
    parser.add_argument("--code", type=int, default=0)
    parser.add_argument("--cap", type=int, default=64)
    parser.add_argument("--declarations-per-part", type=int, default=512)
    parser.add_argument("--name", default="PreparedEJMOriginalAdaptiveCertificate")
    args = parser.parse_args()
    print(json.dumps(generate(Dag(args.dag), args), indent=2))


if __name__ == "__main__":
    main()
