# Fixed census reproduction and source audit

The candidate analyzer reproduces the registered final-layer readings in
both readout orders. The original input is retained unchanged. The audited
version retains eight missing search tracks explicitly and includes one
additional external importer; its measured final counts are also 102/10.
These are finite-table findings, not a registered value function or a
global research policy. No E2 or E3 experiment is included.

## Registration and source binding

[PREREGISTRATION.md](PREREGISTRATION.md) was committed as `6ed88fc1` before
analyzer execution or source auditing. It is an ancestor of the result
commit. The fixed inputs come from experiment commit
`c42809bac1b923b5fc78aaa669125516c3a65ac7`, retained by base
`9567f90e9a8333f0d75d7a04c7e37a57f33dea9c`.

The observation snapshot is main-repository commit
`53ec26edbed43dc4f4862d2db303ad80cb3dbd92`. The analyzer is main-repository
commit `a122da84cd8e7966b67057aef266417aee3ed2df`, built and run in the
designated L3 worktree. That tree's tracked files and branch were unchanged.
The original batch SHA256 is
`8a26793329161a9b0654d5976e187c8e2d167b72cebd273fab57e51d6d90004d`.

## Registered reproduction

| Reading, in each order | Expected | Measured |
| --- | --- | --- |
| Analyzer batch entry exit | 0 | 0 |
| States | 36 | 36 |
| Ordered distinct-state pairs | 1260 | 1260 |
| Final planned-layer inseparable pairs | 102 | 102 |
| Final planned-layer task conflicts | 10 | 10 |
| `residual_sufficiency.sufficient` | false | false |
| `determined_at_layer` | report actual | JSON null |
| Conflict identities | five registered unordered pairs, both orientations | exact match |

The unordered identities are {15354,16702}, {16702,16902}, {16702,17145},
{16702,17605}, {16971,17778}. Null means that no listed layer determines the
task; it does not mean layer zero. The conflict rate is the exact fraction
10/1260.

[analyzer-fixed.json](analyzer-fixed.json) is the actual batch output.
[reconciliation-fixed.json](reconciliation-fixed.json) records 121 passing
comparisons per order. Comparisons include exact pair-set identity and
duplicate checks, not just counts: each unique-capture and leave-one-out
set, each layer's captured set, the final unresolved set, the task-conflict
set, all 1260 first-capture addresses (including nulls), kernel tables,
output-to-kernel reflection, source/plan echoes and exact rates. Pair-set
digests use sorted row-index pairs, whose PR addresses are in the retained
batch. Missing/extra pair lists are empty throughout.

The unchanged `replica.py` reproduces `replica-preview.txt` byte for byte on
the original batch. Its validation and all readings with matching
definitions agree with the analyzer. No registered final-layer
expectation was falsified.

## Catalog scope and planned-readout scope

The analyzer's top-level `escape_pairs`, `unique_capture` and leave-one-out
fields use **all seven catalog roles**, including the task at index 6.
The preview and the plan use **six observation roles**, excluding that
task. This distinction already follows from the candidate's
`Generic.lean` calls to `catalog.fullIndexSet`, `fusedCounts`, and
`selectedAt`; it is not an input correction. The top-level full-catalog
escape count is 92. The planned final unresolved count is 102. The ten
additional pairs are exactly the task conflicts.

| Role | Preview: six-role unique / leave-one-out | Analyzer: seven-role unique / leave-one-out |
| --- | --- | --- |
| shape | 24 / 126 | 14 / 106 |
| admission | 10 / 112 | 6 / 98 |
| resolution | 0 / 102 | 0 / 92 |
| utility | 0 / 102 | 0 / 92 |
| track | 66 / 168 | 44 / 136 |
| period | 12 / 114 | 6 / 98 |
| downstream task | not selected as an observation | 10 / 102 |

The reconciliation checks every analyzer field against its actual role
set and preserves this difference explicitly. Interpreting the top-level
fields as six-role preview fields would be a scope error.

| Layer | Shape-first newly separated / remaining | Period-first newly separated / remaining |
| --- | --- | --- |
| 0 | 0 / 1260 | 0 / 1260 |
| 1 | 616 / 644 | 520 / 740 |
| 2 | 126 / 518 | 336 / 404 |
| 3 | 256 / 262 | 62 / 342 |
| 4 | 0 / 262 (collapsed utility) | 174 / 168 |
| 5 | 148 / 114 | 0 / 168 (collapsed utility) |
| 6 | 12 / 102 | 66 / 102 |

Shape-first adds shape, admission, resolution, utility, track, period.
Period-first adds period, shape, admission, resolution, utility, track.
Every captured pair identity in these layers agrees, as well as the
counts and collapse flags.

## Source audit

[LABEL-AUDIT.md](LABEL-AUDIT.md) displays all seven original and checked
values for each settlement in separate columns with source links.
[labels-audited-v2.json](labels-audited-v2.json) retains their individual
statuses, PR body locators and hashes, immutable source-file locators and
hashes, delivered-module sets and external importer edges.

- 36 rows and 252 fields: 244 source-verified, eight missing, zero checked
  values left inferred. All seven fields are source-verified on 28 rows.
  Verification is of label provenance and repository reads; it is not an
  independent semantic proof-shape, admission or kernel re-verification.
- The eight missing tracks are #11696, #12188, #12833, #13018, #13245,
  #13680, #16361 and #16879. An unlabelled PR and a seat-name suffix do not
  establish quick/deep attribution. Their checked values are JSON null
  with `status=missing`.
- The original inferred shape is #15993. Its PR labels the private helpers
  bind-only and says `result` checks the explicit construction, without
  explicitly classifying that result. The fixed source header at
  `D5/S3/Divergence/QuantumTriangularJensenShannonRefutation.lean:27`
  explicitly says `result: bind-only`. The checked label therefore has a
  direct source, while the original inference remains recorded.
- #14159's original `Proved` was interpreted from its tightness answer.
  Its Scribe explicitly records `ResolutionKind.Proved` at line 83.
- #12833 has six external importers, rather than five. The additional
  source is
  `D5/S3/Quantum/Measurements/BroadcastingChannelExtremaRefutation.lean:309`.
  Its import begins at zero-based character 36829; the original
  `build_census.py` limits each importer read to 30000 characters. The
  boolean task remains true. All 36 boolean task labels agree with the
  source audit.

The audit reads all 8895 D5 Lean blobs at the fixed snapshot, masks nested
comments and string literals, and enumerates 10750 D5 source import edges
without a length cutoff. Added-module identity is determined by each
merge's first-parent diff. Utility is the strongest explicit header among
the added modules. Conclusion types are checked against those modules'
Scribes: 34 explicit resolution claims and two library-result rows.
Admission evidence identifies the new modules, excluding existing owner
re-pins from the aggregate basis. Merge periods are checked against GitHub
merge timestamps.

The importer task is a source import proxy. Live proof-path use is a
separate `not-audited` field throughout and is not substituted for it.

## Audited version and limits

[physics-settlements-36.batch.v2.json](physics-settlements-36.batch.v2.json)
contains the audited labels and importer list. Its output values carry
status alongside value. The eight missing tracks share a missingness
marker; their equality is only equality of missingness, not evidence that
the unknown real tracks agree. The partitions consequently remain the
same in this version, which is a measured result rather than a target.

[analyzer-v2.json](analyzer-v2.json) is the candidate's actual output for
this version. [reconciliation-v2.json](reconciliation-v2.json) records 121
passing comparisons per order: 36 states, 1260 ordered pairs, 102 planned
unresolved pairs, ten conflicts with the same identities, false
sufficiency and null first determining layer. The original preview's
string-only fiber formatter cannot display null track values; v2 is
checked directly against its table presentation, preserving nulls.

Neither the analyzer nor this audit proves that source labels are correct
semantic judgments. The analyzer retains its unverified-source boundary.
The shape/exposure confounding and the distinction between import and
proof use remain. No value function is registered, including no claim
that an unregistered function is a zero vector.

## Reproduction and validation

From the experiment repository root, with Python 3.12 available:

```sh
python3.12 -B docs/reports/physics-settlement-value-census/audit_census.py \
  --main-repo /Users/auricstudio/trureturing-0addbbde-l3 \
  --out docs/reports/physics-settlement-value-census
python3.12 -B docs/reports/physics-settlement-value-census/reconcile.py \
  docs/reports/physics-settlement-value-census/physics-settlements-36.batch.json \
  docs/reports/physics-settlement-value-census/analyzer-fixed.json \
  --registered --out docs/reports/physics-settlement-value-census/reconciliation-fixed.json
python3.12 -B docs/reports/physics-settlement-value-census/reconcile.py \
  docs/reports/physics-settlement-value-census/physics-settlements-36.batch.v2.json \
  docs/reports/physics-settlement-value-census/analyzer-v2.json \
  --out docs/reports/physics-settlement-value-census/reconciliation-v2.json
```

The source audit obtains PR metadata using read-only `gh pr view`; edited
body locators fail when their selected labels no longer match. To rerun
the analyzer, use the registered build and batch entry in the designated
L3 tree and redirect the relevant batch from this directory. Select Python
3.12 in that command's PATH; the host's default Python lacks `tomllib`.
Do not replace the fixed original batch with v2.

| Check | Exit |
| --- | --- |
| `make -C tools auric-fib-analysis-build` in the designated L3 tree | 0 |
| Actual `tools/scripts/auric-fib-analysis.sh --batch`, original batch | 0 |
| Actual `tools/scripts/auric-fib-analysis.sh --batch`, audited batch | 0 |
| `audit_census.py`, fixed git snapshot and all 36 PRs | 0 |
| `reconcile.py --registered`, original batch and stored preview | 0 |
| `reconcile.py`, audited batch | 0 |

Only the necessary programs, data, source evidence and findings are
published. Process logs, raw PR response copies and runner receipts are
excluded. This delivery claims no separate reviewer consensus.
