# Bounded chordal graph polytope probe

This program searches for an exact finite counterexample to the full original
Linusson–Restadh Conjecture 3.4. A result without a candidate leaves that
conjecture unresolved. Computed candidates require separate independent
checking and subsequent Lean verification; this program supplies neither.

Run from the experiments repository root with the same interpreter:

```sh
/usr/bin/python3 experiments/rhombus-conj3.4-23531-ownedpath/probe.py experiments/rhombus-conj3.4-23531-ownedpath/registration.json
```

The original attempt is single use. Its `execution-start.json` prevents a second
execution in this directory. Do not remove it or restart/resume this attempt.
The command above documents the original invocation. Later independent checking
must use the saved data under separate authorization and its own bounds.

Dependencies: Python 3.9 or later, standard library only, on a POSIX host with
`signal.setitimer`. Actual interpreter and platform are in `results.json`.
There are no required third-party Python packages or solver installations.

`registration.json` retains the functional historical input: one complete,
unchanged REST-verified preregistration comment body (including the exact event
and both predictions), its immutable provenance, the operator timing and the
fields consumed by the unchanged program. `strict_preregistration` retains only
the consumed deadline; the complete prediction block has one source of truth in
`comment.body`. Unconsumed transport/gate receipts have been removed.
`execution-start.json` is now a minimal persistent used marker; its existence
continues to prohibit another execution. Original timing and caps remain in
`results.json`. The probe makes no network or Git calls.

The fixed route uses n=7, at most 65,536 raw proposals, at most 256 completely
processed eligible pairs, 3 through 12 Delta components, all 2^k selections
(up to 4,096) for each pair, and 3,600 seconds of actual computation. Preparation
has a separate 1,200-second allowance within the 4,800-second overall allowance.
The first applicable cap stops the attempt, with no reseeding or replacement.
All graph and solver validation occurs within the actual computation timer.

For decimal j starting at zero, the proposal is SHA256 of the exact prefix
`Linusson-Restadh/2305.05275v1/Conj3.4/exact-refutation/v1/` immediately followed
by j, with no whitespace or leading zeros. The first six digest bytes are an
unsigned big-endian integer reduced modulo 2^42. Low/high 21-bit edge masks are
ordered numerically; equal and duplicate pairs are counted and discarded.
Edge bit positions follow lexicographic unordered pairs of labels 0 through 6.
No quotient by relabelling is taken. Coordinate bit positions follow increasing
numeric seven-bit subset masks with at least two set bits (120 coordinates).

Chordality is checked by the complete finite list of injective cycles of lengths
4 through 7. Fixing the smallest vertex first and the orientation by second <
last removes only repeated descriptions of the same cycle. A cycle obstruction
requires exactly the cyclic edges in its induced subgraph, excluding every
chord. Known-case exclusions are nested edge sets, two split graphs (each tested
by all clique/independent-set partitions), and a changed vertex neighborhood
whose two neighborhoods are both complete in both graphs. The latter checks all
four completeness conditions of source Lemma 3.15.

Delta has every differing clique subset as a vertex. Two distinct vertices are
adjacent exactly when their intersection is another differing clique subset.
All its connected components are retained. Each selection preserves every
common coordinate and assigns one value to each whole component. The unique
graph determined by its size-two coordinates is checked for both chordality and
exact equality of every clique indicator. Every processed pair file retains the
complete compatible family and a concrete obstruction for each failed selection.

A new complementary compatible pair is verified in all 120 coordinates and
excludes both original endpoints. For pairs without one, rational phase I
simplex with Bland's rule seeks a convex combination of the other actual
compatible graphs on a strictly interior point of the endpoint segment. Solver
validation uses three finite linear systems within the same execution timer;
validation and pair LP calls are counted separately. No floating arithmetic is
used for mathematical coordinates or weights. Any returned hit is checked again
for actual graph chordality, endpoint exclusions, rational strict interiority,
normalization, nonnegativity and equality in all 120 coordinates. Negative
solver results are search outcomes and are not exported as face certificates.

`results.json` contains the final counters, cap, timing, dependency information
and candidate if present. `pairs/pair-J.json` contains exact native edge masks,
components, complete families, invalid-selection cycle/clique witnesses and any
rational certificate. Pair files are individually hashed in `results.json`.
A current partial pair, if interrupted, is explicitly distinguished from a
completely processed pair. Atomic writes preserve the latest recoverable data.
No child processes are launched by the probe.

The remaining formal obligations are the source/native graph and full imset
correspondence, actual real convex-hull vertex correspondence, cycle enumeration
completeness, complementary-component witness completeness, and the ordinary
face consequence of an exact interior convex equality. Neither program bytes,
solver statuses nor a finite absence of hits discharge those obligations.

## Completed result

The single fixed n=7 probe stopped at the first cap of 256 fully processed eligible pairs (6,989 raw proposals, last index 6,988; 7,480 complete component selections). It retained 120 pairs with new rhombus witnesses and 136 without a new complementary pair; 136 exact rational phase-I LP calls returned no feasible point. No exact counterexample certificate was found. These negative LP outcomes are not face certificates. The original all-positive-n Conjecture 3.4 and the all-pairs n=7 claim remain unresolved.

Execution started at 2026-10-29T06:25:12.820185+00:00 and ended at
2026-10-29T06:25:16.831085+00:00; recorded elapsed time was 4.010819625 seconds.
No second attempt, proposal stream, or LP search was performed after the cap.
The no-hit supplies no general implication and qualified no successor candidate.
Resolution KPI credit is zero. Native graph/imset, real-hull vertex, finite-cycle,
component/witness and ordinary-face bridges remain uncompiled.

A separate completed evaluator assessed the bounded experiment as
valid-unresolved/defer and approved factual unresolved reporting. The
`formal_status` inside unchanged `results.json` is the original probe's
historical output, captured before that evaluation; no Lean proof was produced.
The original program, all 256 pair files, results, source references and licenses
are preserved byte-for-byte. `manifest.json` indexes the retained files, with
byte counts and SHA256 hashes; the manifest does not hash itself.
