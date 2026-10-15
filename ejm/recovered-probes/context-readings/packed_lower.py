"""Prepare exact paged scalar lower lookups; never edit live Lean inputs.

Literal values are untrusted until the complete finite kernel equalities
against every original mirror lower succeed. No new global acceptance claim.
"""

from pathlib import Path
import re
import struct

ROOT = Path(__file__).resolve().parents[1]
DATA = Path("/Users/macstudio/omega/trureturing-01a12631-certificate/.ejm-proof-data/rounded32")
OUT = Path(__file__).resolve().parent / "packed-lower-candidates"
OUT.mkdir(exist_ok=True)
NS = "D5.S3.Quantum.Entanglement.EJM"
SCALE = 6987435
PAGE = 128
BIAS = 1 << 63
WORD = 1 << 64


def rows(block, count, dimension):
    raw = (DATA / f"local-dag-{block}.bin").read_bytes()
    assert struct.unpack_from("<I", raw)[0] == count
    off = 4
    values = []
    for _ in range(count):
        _, lower, split = struct.unpack_from("<IqB", raw, off)
        off += 13 + (16 if split < dimension else 0)
        scaled = lower * SCALE
        assert -BIAS <= scaled < BIAS
        values.append(scaled)
    assert off == len(raw)
    return values


def pages(values, stem):
    words = [sum((value + BIAS) << (64 * k)
                 for k, value in enumerate(values[start:start + PAGE]))
             for start in range(0, len(values), PAGE)]
    # Check the encoding completely as a producer diagnostic. Only the Lean
    # receipt below earns equality with the original source.
    assert [((words[j // PAGE] >> (64 * (j % PAGE))) % WORD) - BIAS
            for j in range(len(values))] == values
    parts = []

    def tree(lo, hi):
        if lo == hi:
            return ".empty"
        mid = (lo + hi) // 2
        return f"(.node {mid} 0x{words[mid]:x} {tree(lo, mid)} {tree(mid + 1, hi)})"

    def split(lo, hi):
        if lo == hi:
            return ".empty"
        mid = (lo + hi) // 2
        if hi - lo <= 32:
            term = tree(lo, hi)
        else:
            term = f"(.node {mid} 0x{words[mid]:x} {split(lo, mid)} {split(mid + 1, hi)})"
        name = f"{stem}Part{len(parts)}"
        parts.append(f"@[noinline] private def {name} : RawTree :=\n  {term}\n")
        return name

    root = split(0, len(words))
    return "\n".join(parts) + f"\nprivate def {stem} : RawTree := {root}\n"


def cases(names):
    heads = [f"Fin.cases {name} (" for name in names]
    lines = ["  " + "".join(heads[i:i + 16]) for i in range(0, len(heads), 16)]
    lines.append("  (fun b : Fin 0 => Fin.elim0 b)" + ")" * len(heads))
    return "\n".join(lines)


values = {"triangle": rows(0, 41294, 10), "singleton": rows(1, 8689, 9)}
data_source = f"""/- Prepared uncompiled scalar representation candidate.
Every signed 64-bit value is an untrusted literal. Full finite equality with
every original MirrorRow.lower is required before use in global replay.
128 original rows per packed page; the final short pages are retained.
No state, transition, score, gap, assignment or value branch is changed. -/

import {NS}.EJMPartialBoundData

set_option autoImplicit false
set_option relaxedAutoImplicit false
set_option maxRecDepth 100000
set_option linter.style.longLine false

namespace {NS}.EJMPackedLowerData
open EJMPartialBoundData (RawTree)

"""
for stem, count in [("triangle", 41294), ("singleton", 8689)]:
    data_source += pages(values[stem], stem + "Pages")
    data_source += f"""
def {stem}Lower (j : Fin {count}) : ℤ :=
  (({stem}Pages.get (j.val / {PAGE}) /
    2 ^ (64 * (j.val % {PAGE})) % 18446744073709551616 : ℕ) : ℤ) - 9223372036854775808

"""
data_source += f"end {NS}.EJMPackedLowerData\n"
(OUT / "EJMPackedLowerData.lean").write_text(data_source)

receipt_source = f"""/- Prepared uncompiled complete finite scalar receipts.
Every actual row is compared with the original complete mirror's lower field.
The two domain totals are 41,294 and 8,689; no producer minimum is a premise. -/

import {NS}.EJMPackedLowerData
import {NS}.EJMTriangleCacheData
import {NS}.EJMSingletonCacheData

set_option autoImplicit false
set_option relaxedAutoImplicit false
set_option Elab.async false
set_option maxRecDepth 100000
set_option maxHeartbeats 100000000
set_option trace.profiler true
set_option linter.style.longLine false

namespace {NS}.EJMPackedLowerReceipts
open EJMPackedLowerData

"""
for stem, count, title in [("triangle", 41294, "Triangle"), ("singleton", 8689, "Singleton")]:
    names = []
    for start in range(0, count, PAGE):
        length = min(PAGE, count - start)
        index = f"⟨{start} + k.val, by have hk := k.isLt; omega⟩"
        name = f"{stem}_{start}"
        names.append(name)
        receipt_source += (f"private theorem {name} : ∀ k : Fin {length}, "
                           f"{stem}Lower {index} = "
                           f"(EJM{title}CacheData.mirror {index}).lower := by decide +kernel\n")
    total = len(names)
    receipt_source += f"""
private def {stem}Length (b : Fin {total}) : ℕ := min {PAGE} ({count} - {PAGE} * b.val)
private def {stem}Index (b : Fin {total}) (k : Fin ({stem}Length b)) : Fin {count} :=
  ⟨{PAGE} * b.val + k.val, by
    have hb := b.isLt
    have hk := k.isLt
    dsimp [{stem}Length] at hk
    omega⟩

private theorem {stem}_blocks : ∀ b : Fin {total}, ∀ k : Fin ({stem}Length b),
    {stem}Lower ({stem}Index b k) = (EJM{title}CacheData.mirror ({stem}Index b k)).lower :=
{cases(names)}

theorem {stem}_eq (j : Fin {count}) :
    {stem}Lower j = (EJM{title}CacheData.mirror j).lower := by
  let b : Fin {total} := ⟨j.val / {PAGE}, by have hj := j.isLt; omega⟩
  let k : Fin ({stem}Length b) := ⟨j.val % {PAGE}, by
    have hj := j.isLt
    dsimp [{stem}Length, b]
    omega⟩
  have he : {stem}Index b k = j := by
    apply Fin.ext
    dsimp [{stem}Index, b, k]
    omega
  simpa only [he] using {stem}_blocks b k

#print axioms {stem}_eq

"""
receipt_source += f"end {NS}.EJMPackedLowerReceipts\n"
assert len(receipt_source.splitlines()) < 800
(OUT / "EJMPackedLowerReceipts.lean").write_text(receipt_source)

bridge = f"""/- Prepared uncompiled exact scalar aggregate bridge.
The full original cache states and future transitions remain unchanged.
Only delta's next lower projection uses the fully checked scalar table;
the current row is already needed for its unchanged transition payload.
There is no replay acceptance theorem in this candidate. -/

import {NS}.EJMPackedLowerReceipts
import {NS}.EJMAggregateCachePacking

set_option autoImplicit false
set_option relaxedAutoImplicit false
set_option Elab.async false

namespace {NS}.EJMPackedLowerAggregate
open EJMCachePacking
open EJMLocalCacheData
open EJMPackedLowerData
open EJMAggregateCachePacking (AggregateCache)

def triangleDelta (g : Fin 20) (a : Fin 4) (k : Fin 16) (word : ℕ) : ℤ :=
  match triangleIncidence g k with
  | none => 0
  | some i =>
    let state := decodeState (by decide : 0 < 41294) word
    let current := EJMTriangleCacheData.mirror state.1
    let next := readTransition (by decide) current (trianglePreimage state.2 i).val
      (EJMPartialBoundData.outcome state.2.2 a)
    triangleLower next.1 - current.lower

def singletonDelta (g : Fin 20) (a : Fin 4) (k : Fin 16) (word : ℕ) : ℤ :=
  match singletonIncidence g k with
  | none => 0
  | some i =>
    let state := decodeState (by decide : 0 < 8689) word
    let current := EJMSingletonCacheData.mirror state.1
    let next := readTransition (by decide) current (singletonPreimage state.2 i).val
      (EJMPartialBoundData.outcome state.2.2 a)
    singletonLower next.1 - current.lower

theorem triangleDelta_eq (g : Fin 20) (a : Fin 4) (k : Fin 16) (word : ℕ) :
    triangleDelta g a k word = EJMAggregateCachePacking.triangleDelta g a k word := by
  cases hi : triangleIncidence g k with
  | none => simp only [triangleDelta, EJMAggregateCachePacking.triangleDelta, hi]
  | some i =>
    simp only [triangleDelta, EJMAggregateCachePacking.triangleDelta, hi,
      EJMPackedLowerReceipts.triangle_eq, EJMCacheState.triangleStep]

theorem singletonDelta_eq (g : Fin 20) (a : Fin 4) (k : Fin 16) (word : ℕ) :
    singletonDelta g a k word = EJMAggregateCachePacking.singletonDelta g a k word := by
  cases hi : singletonIncidence g k with
  | none => simp only [singletonDelta, EJMAggregateCachePacking.singletonDelta, hi]
  | some i =>
    simp only [singletonDelta, EJMAggregateCachePacking.singletonDelta, hi,
      EJMPackedLowerReceipts.singleton_eq, EJMCacheState.singletonStep]

def delta (cache : PackedCache) (g : Fin 20) (a : Fin 4) : ℤ :=
  EJMStreamingCachePacking.sumWords 16 (triangleDelta g a) cache.triangle +
    EJMStreamingCachePacking.sumWords 16 (singletonDelta g a) cache.singleton

theorem delta_eq (cache : PackedCache) (g : Fin 20) (a : Fin 4) :
    delta cache g a = EJMAggregateCachePacking.delta cache g a := by
  unfold delta EJMAggregateCachePacking.delta
  have ht : triangleDelta g a = EJMAggregateCachePacking.triangleDelta g a := by
    funext k word
    exact triangleDelta_eq g a k word
  have hs : singletonDelta g a = EJMAggregateCachePacking.singletonDelta g a := by
    funext k word
    exact singletonDelta_eq g a k word
  rw [ht, hs]

def step (cache : AggregateCache) (g : Fin 20) (a : Fin 4) : AggregateCache :=
  ⟨EJMStreamingCachePacking.step cache.words g a,
    cache.total + delta cache.words g a⟩

theorem step_eq (cache : AggregateCache) (g : Fin 20) (a : Fin 4) :
    step cache g a = EJMAggregateCachePacking.step cache g a := by
  simp only [step, EJMAggregateCachePacking.step, delta_eq]

#print axioms triangleDelta_eq
#print axioms singletonDelta_eq
#print axioms delta_eq
#print axioms step_eq

end {NS}.EJMPackedLowerAggregate
"""
(OUT / "EJMPackedLowerAggregate.lean").write_text(bridge)
for stem, xs in values.items():
    print(stem, "signed64-safe rows", len(xs), "min", min(xs), "max", max(xs))
for path in sorted(OUT.glob("*.lean")):
    print(path.relative_to(ROOT), "lines", len(path.read_text().splitlines()),
          "bytes", path.stat().st_size)
