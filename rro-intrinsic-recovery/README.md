# Intrinsic recovery from the actual paired-power sources

`exact_checks.py` is a standalone exact finite program; `exact-results.json`
is its deterministic result. It uses Python 3.12+ and the standard library,
with Git access to the pinned helper objects. Assertions must remain enabled.
There are no third-party packages, native ports, new source draws or main-repo
builds. Validation covers the stated finite cases through the author's exact run and
alternate-hashseed comparison; it does not constitute independent review or
Lean certification.

The mathematical input is trureturing revision
`d6e6417f0e055ac0701e5853b9eaf70a316ef709`: UJB Proposition 3.2 and
§§5,9,12,15–18, POSITIVE_RELATION_CERTIFICATES §§1–7, and
FINITE_GROUP_POWER_MAP_SEPARATION §§1–3. The experiment base is
`b490f4f17c97645ce7c6495ca8635bce3d3a1af7`. The program loads the existing
literal helper at `865f38257ead1b4e857ff915c23b0b04276bd7f3` and the positive
certificate/tree/read helpers at the experiment base, in memory. Their exact
paths and SHA256 hashes are in the result. Their unrelated checks are not run.

Run from this worktree, using the required cooperative entry:

```sh
python3.12 -B /Users/auricstudio/trureturing/tools/scripts/worktree/worktree_protocol.py \
  --source /Users/auricstudio/trureturing-experiments with --path "$PWD" \
  --read rro-structural-recovery --read rro-positive-certificates \
  --write rro-intrinsic-recovery -- \
  python3.12 -B rro-intrinsic-recovery/exact_checks.py \
  --output rro-intrinsic-recovery/exact-results.json
```

For a read-only deterministic comparison, replace the write scope with
`--read rro-intrinsic-recovery` and replace `--output` with `--check`.
The JSON omits timing and interpreter patch versions so identical mathematical
inputs give identical bytes. Arithmetic uses integers, F3 and exact fractions.
Labels are MSB-first binary integers in 7/8 bits. Program indices are zero-based;
reported native sample indices explicitly marked `1based` are one-based.

The program reconstructs all 4,374 elements of each actual group S_i, checks
inversion, classifies orders, and obtains the 1,944 literal edges from squares
and cubes of the order-six elements. It compares this adapter with the pinned
literal relation, including all 72/243 endpoints, all w3 and all nine z values.
The graphs G_i themselves are not groups.

`Geometry` receives only certified transported literal edges and their labels.
It computes equal-neighborhood classes. For G0 these are 36 pairs on the left
and 27 nine-element classes on the right, with 108 quotient incidences.
Disjoint quotient lines have 0,2,3 common transversals in respectively
36,324,108 cases. Connecting pairs with three transversals recovers four
parallel classes of nine lines, and the transversal unions give 18 marked
nine-point planes. Intersections construct affine charts from an origin,
three ordered marked directions and one of two orientations on each axis.
The resulting origin-free operation is H(p,q,r)=p-q+r.

Every one of the 5,184 chart choices is checked against the original source
coordinates under each of four declared full relabelings: identity, reversed
ambient labels, two odd affine permutations of ambient labels, and SHA256
ordering. The full maps are recorded. There are 20,736 charts and 559,872
chart-point comparisons. For each relabeling all 19,683 H triples are compared
directly. All charts are verified affine point by point; cancellation in
(b+Au)-(b+Av)+(b+Aw) then proves H agrees across their triples, including
choices not used by the canonical implementation. Its full incidence has
117 lines, 13 parallel classes and 39 planes; point/line/plane incidences are
checked. The four marked directions and actual endpoint-label associations
lift back to exactly E, checked on all 17,496 literal endpoint pairs.

These constructions commute with relabeling: equal neighborhoods, intersection,
transversal counts and connected components are preserved by bijections.
The finite choice coverage therefore transports to arbitrary full injections;
the computation enumerates only the four stated injections. No full-injection
exhaustiveness, general graph-recognition or new rigidity theorem is claimed.
The Cara–Rottey–Van de Voorde rigidity supplier is not needed for this finite
construction, and no literature novelty or Lean certification is asserted.

Actual group inversion sends an edge (a,t) to (a^-1,t). G0's unique other left
twin supplies a^-1's label, so certified E and an acquired actual input edge
recover its inverse labels. Outside that input domain the program refuses.
G1 has nine K8,27 components, retaining all free (w3,z) choices, and is refused
by the G0 constructor. For each relabeling the program checks 432 G1 witnesses:
swap an input's inverse partner with any of six other twins while fixing the
input. Change the full map once and transport every hidden prepared edge.
The same E and all 1,945 prepared label pairs remain identical, while that
same input's inverse label changes. These are positive realizations whenever
G1 has positive prior weight. For the odd-permutation example the acquired
input is labels (11,19), but its inverse labels are (88,19) or (37,19).

The unresolved G1 reference is a shared matching. All 40,320 local bijections
on one eight-element class yield 105 matchings, each 384 times. Given E and
any supported prepared label array, a specified input's partner has mass 1/7
on each other twin. With only its left label acquired, inverse edges have
189 possibilities of mass 1/189. Mutual partners have joint probability 1/7;
independent reference redraws would give 1/49. Repeated targets agree surely.
Separately, inversion pushes an unread uniform edge to a uniform edge in both
sources. That marginal statement supplies neither pointwise G1 recovery nor
permission to draw an inverse label or redraw the latent matching.

On the certified-E domain the program reuses PCR.9/UJB's indexed law:
each original slot is uniform on its own E intersected with its acquired bit
rectangle, and different slots are conditionally independent. G0 inverse
pushforwards agree with the actual group inverse for 1,979 positive masks per
relabeling: empty, all 30 single-bit rectangles, four representative field/mixed
rectangles, and all 1,944 complete edge rectangles. Sequential partial
completion telescopes to 1/27, including a repeat with probability one.
H output laws count literal lifts per point, rather than replacing them by
independent quotient samples. The three-unread case counts 1,944^3 literal
arrays via all 19,683 point triples; partial/complete/unread and repeated-index
cases are also exact. Reusing one partial index gives point agreement one;
an independent redraw would give 1/3. Same-index inverse agreement is one,
whereas two distinct equally restricted slots agree with probability 1/27.

One common public controller scans the first 576 slots of a fixed N=580 native
preparation in either source. The G1 fixture pads its 324-pair certificate by
repeated witnessed edges; these slots and draws still count. Every acquisition
uses the existing stationary original-root Read interpreter. Both scans cost
8,640 invocations, operand-length sum 37,337,760, maximum length 8,641.
Original slot 577 is partial, 578 complete, and 579–580 remain unread. The result
retains every original indexed cache as compact [mask,value] integer rows;
integer bit h encodes local coordinate h, with no merging of repeated samples. In G0
the controller completes slot 577 and computes inverse labels (31,2) and an
H point represented by its nine actual right labels. G1 reports no pointwise
inverse output. Full costs are respectively 8,674/8,666 invocations and
8,672/8,664 distinct addresses; both have maximum operand length 8,671.
Preparation costs 580 paid draws and 43,504 native events in either fixture.
These are supported realizations of the original random preparation, not a
source-aware acquisition policy or an acquisition probability improvement.

Invocation count charges repeats/constants/absences. The distinct-address
tariff of ActualTreeReadoutAcquisition is reported separately. Literal operand
length, cache/label storage, incidence reconstruction, finite arithmetic,
controller state, preparation, permissions and external randomness remain
separate resources; none is supplied free by the output law. Law realization
means acquiring missing bits of original slots and applying deterministic
postprocessing. It is not an additional group-operation or sampling port.

Twenty-six nonzero translations preserve G0 observations but change the hidden
origin, and a concrete right-twin swap changes z without changing E. Swapping
two completed old slots preserves their compatible-generator likelihood but
changes indexed replay. NativePromiseCausalBoundary supplies the independent
serial-history limitation: opposite orders of the legal insertions at two
disjoint inactive C slots reach the same endpoint. Thus the endpoint can retain
its labelled generation partial order and grade without retaining the selected
serial order. InvolutionDescent, TargetRecoveryCriterion, PushforwardComposition
and RRO Proposition 92.2 are inspected suppliers for descent, fiber constancy
and law composition, not fresh compiled applications here.

Only the specified chart choices, literal source domains, local bijections and
stated finite families are exhaustive. Mixed masks, relabelings and native
controllers/words are representative. Arbitrary indexed continuation and
positive certificate soundness use the existing suppliers. Closure inverses
retain their fixed public prior, fixed indexed partial cache and at least one
wholly unread original sample; the unread empty-product law has no general
inverse. G1 component permutations preserve incidence but need not preserve
any chosen affine-plane heap, as all 36 component transpositions demonstrate.
No hidden injection/origin, actual serial history, source mutation, native
random/quantum operation, optimality or resource improvement is recovered.
