"""Prepare compact 128-row receipt sources without touching live Lean inputs.

Every existing finite kernel proposition and its decide proof is retained.
Dependent Fin.cases dispatches exactly one existing chunk for each block;
no arithmetic conclusion is transferred across different row contexts.
The generated candidates have not been compiled or accepted.
"""

from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = Path(__file__).resolve().parent / "compact-candidates"
OUTPUT.mkdir(exist_ok=True)


def cases(names):
    # Sixteen syntactic Fin.cases heads per physical line. The final empty
    # domain is eliminated, so no finite case is silently omitted.
    pieces = [f"Fin.cases {name} (" for name in names]
    lines = ["  " + "".join(pieces[i:i + 16])
             for i in range(0, len(pieces), 16)]
    lines.append("  (fun b : Fin 0 => Fin.elim0 b)" + ")" * len(names))
    return "\n".join(lines)


for kind, count in [("Triangle", 41294), ("Singleton", 8689)]:
    name = f"EJM{kind}CacheReceipts"
    original = ROOT / "D5/S3/Quantum/Entanglement/EJM" / f"{name}.lean"
    source = original.read_text()
    chunks = re.findall(
        r"private theorem ((?:mirror|transitions)_[0-9]+)\s*:.*?"
        r":= by\s+decide \+kernel", source, re.S)
    raw_chunks = re.findall(
        r"private theorem (?:mirror|transitions)_[0-9]+\s*:.*?"
        r":= by\s+decide \+kernel", source, re.S)
    block_count = (count + 127) // 128
    assert len(chunks) == len(raw_chunks) == 2 * block_count
    mirror_names = [x for x in chunks if x.startswith("mirror_")]
    transition_names = [x for x in chunks if x.startswith("transitions_")]
    starts = [128 * b for b in range(block_count)]
    assert mirror_names == [f"mirror_{i}" for i in starts]
    assert transition_names == [f"transitions_{i}" for i in starts]
    header = source[:source.index("set_option maxHeartbeats 100000000 in")]
    header = re.sub(r"Prepared source: not yet compiled\. ", "", header)
    header += "set_option maxHeartbeats 100000000\n"
    header += "set_option trace.profiler true\n"
    header += "set_option linter.style.longLine false\n\n"
    compact = header + "\n".join(re.sub(r"\s+", " ", x) for x in raw_chunks) + "\n\n"
    compact += f"""private def chunkLength (b : Fin {block_count}) : ℕ :=
  min 128 ({count} - 128 * b.val)

private def chunkIndex (b : Fin {block_count}) (k : Fin (chunkLength b)) : Fin {count} :=
  ⟨128 * b.val + k.val, by
    have hb := b.isLt
    have hk := k.isLt
    dsimp [chunkLength] at hk
    omega⟩

private theorem mirror_blocks : ∀ b : Fin {block_count}, ∀ k : Fin (chunkLength b),
    checkMirrorFields (chunkIndex b k) = true :=
{cases(mirror_names)}

private theorem transitions_blocks : ∀ b : Fin {block_count}, ∀ k : Fin (chunkLength b),
    checkTransitions (chunkIndex b k) = true :=
{cases(transition_names)}

private theorem mirror_fields_accepted (j : Fin {count}) : checkMirrorFields j = true := by
  let b : Fin {block_count} := ⟨j.val / 128, by have hj := j.isLt; omega⟩
  let k : Fin (chunkLength b) := ⟨j.val % 128, by
    have hj := j.isLt
    dsimp [chunkLength, b]
    omega⟩
  have he : chunkIndex b k = j := by
    apply Fin.ext
    dsimp [chunkIndex, b, k]
    omega
  simpa only [he] using mirror_blocks b k

"""
    mirror_bridge = source[source.index("theorem mirror_accepted "):]
    mirror_bridge = mirror_bridge[:mirror_bridge.index("\nset_option maxHeartbeats")]
    compact += mirror_bridge + "\n\n"
    compact += f"""theorem transitions_accepted (j : Fin {count}) : checkTransitions j = true := by
  let b : Fin {block_count} := ⟨j.val / 128, by have hj := j.isLt; omega⟩
  let k : Fin (chunkLength b) := ⟨j.val % 128, by
    have hj := j.isLt
    dsimp [chunkLength, b]
    omega⟩
  have he : chunkIndex b k = j := by
    apply Fin.ext
    dsimp [chunkIndex, b, k]
    omega
  simpa only [he] using transitions_blocks b k

"""
    compact += source[source.index("theorem assignment_eq "):]
    # Static source checks only. Actual Lean elaboration and kernel validation
    # remain mandatory after the running complete receipt reaches a terminal.
    compact_chunks = re.findall(
        r"private theorem (?:mirror|transitions)_[0-9]+\s*:.*?"
        r":= by\s+decide \+kernel", compact, re.S)
    assert [re.sub(r"\s+", " ", x) for x in raw_chunks] == compact_chunks
    assert len(compact.splitlines()) <= 1000
    out = OUTPUT / f"{name}.lean"
    out.write_text(compact)
    print(out.relative_to(ROOT), "lines", len(compact.splitlines()),
          "exact original chunks", len(raw_chunks), "domain", count)
