# C1 finite packing evidence

C1 remains unresolved. This directory retains reusable finite packing code,
existing data, and an independent replay of the saved evidence. It supplies no
whole-source proof or counterexample, no supported successor solving route, and
no Lean/kernel closure, admission, completion, novelty claim, or KPI.

## Source and boundary

Alexandr Kostochka et al., *Turán number of three disjoint cliques of arbitrary
sizes*, [arXiv:2610.38430v1, Section 11, C1](https://arxiv.org/html/2610.38430v1#ThmConjecture1),
2026. Source license: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
The source attribution and license are supplied documentary evidence; this
retention replay did not fetch the article or repeat literature/library checks.
The source statement is paraphrased below and in `encoding.md`; the article is
not redistributed here. See `NOTICE.md` for licensing and modification details.

For the fixed conventional event, for every integer ell >= 2, sizes with
each earlier k_i >= k_ell >= 2, n >= sum(k), and simple undirected n-vertex graph
G, put s = n - sum(k) + 1. If G has disjoint independent sets of the first
ell - 1 exact sizes and e(G) < (2ell - 1)s, C1 asks for disjoint independent
sets of all ell exact sizes. Every set may be chosen again. Edges between sets
are permitted. The earlier sizes have no ordering requirement among themselves.

C1's source block does not explicitly give a lower bound on ell. The possible
ell = 1 convention remains unresolved for a positive whole-source conclusion;
C2's explicit ell >= 4 restriction cannot be imported into C1. Here ell >= 4
is solely the retained finite domain. A valid refutation in that domain would
be within both readings, but none is present.

## Validated extent

`validation-results.json` records deterministic results from `verify_retained.py`:

- All 30,000 saved inputs have valid full-packing witnesses in the original
  size order. A separate endpoint/set checker verifies graph endpoints,
  symmetry, irreflexivity, strict source edge bound, exact sizes, disjointness,
  and independence. Their first ell - 1 sets also witness the partial premise.
  Inputs may repeat graphs; 30,000 counts input records, not distinct graphs.
- There are 7,500 records for each seed 261038430, 1, 2, and 3. The finite
  domain is ell in {4,5,6}, k_ell in {2,3}, each earlier size in
  {k_ell,k_ell+1,k_ell+2}, t in {0,...,6}, and n = sum(k) + t.
- The original small controls replay with 4,096 literal-versus-naive and
  4,096 capacity-coloring-versus-naive comparisons on every labeled n=5
  graph and four size vectors; 1,751 clique-versus-literal comparisons on the
  original n<=10 partition domain; fixed-partial failure with full rechoice
  success; and the K7 plus three isolates equality obstruction. Equality is
  outside the strict hypothesis and is not a C1 counterexample.
- The saved clique table has 644 parameter rows. Its stored counts sum to
  778,206 parameter/partition pairs and 344,943 strict-eligible pairs. Integer
  partition counts are independently checked by dynamic programming. The
  644 saved obstruction records have correct component/edge/bound arithmetic
  and partial/full cut classifications. The large sweep, all its strict-edge
  classifications, and its claims of minimum edge counts were not replayed.

The fixed input stream was reconstructed once from the retained deterministic
generator to make the existing data explicit; this took no packing search and
matched the previously supplied stream hash. Default verification reads this
fixed stream and does not call the generator. The direct witness checker does
not call either packing solver on the 30,000 inputs. No coverage was extended.

The independently reproduced canonical input stream SHA-256 is
`aebdf7d53ba8177a43a767ec3d4809156fe2ed083c52748b5f28e57696aab875`.
The witness stream SHA-256, in input-stream order, is
`3c45d8afd07cc5748a2d694c935b646aff650751151177240a73c521bfec93dd`.
These are canonical JSON-line stream digests, distinct from compressed-file and
saved JSON-file digests in `SHA256SUMS`.

## Replay

Run from this directory with Python 3.13; only the standard library is needed:

```sh
python3.13 -B verify_retained.py --check-report validation-results.json
shasum -a 256 -c SHA256SUMS
```

The checker raises on any failed obligation. `NodeLimit` means unknown, never
nonpacking. The solvers are used only for the original small controls. The
capacity checker retains the existing branching/cache algorithm with its
process and wall-time wrappers replaced by an explicit node bound.

## Files and limits

`packing_core.py`, `packing-witnesses.json`, and `clique-results.json` retain
their original bytes. `upstream-artifact-sha256.json` records the necessary
original artifact fingerprints, including the sources of the cleaned encoding
and extracted capacity function. `fixed-inputs.jsonl.gz` contains only the
existing input records in their original generation order, with sorted edge
lists and canonical JSON lines. `capacity_packing.py`, `verify_retained.py`,
`encoding.md`, and `validation-results.json` make the replay portable.

The whole-target gap remains a proof for arbitrary source-admissible graphs,
sizes, n and ell, or a finite source-valid graph with complete exclusion of
every full packing allowing all sets to change. A clique-family allocation
test supplies no reduction from arbitrary graphs. Finite positive witnesses
cannot establish the universal event. The result remains unresolved, neither
false nor scored. All six assessment coordinates remain separate: generativity,
new connection gain, general sharpness, and whole-target generality gain are
unmeasured; discrimination is not re-evaluated; continuation has no supported
successor. Unmeasured does not mean zero, and these are no conjunctive gate.

No process plans, progress records, forecasts, dialogue, logs, status receipts,
or repository snapshots are retained.
