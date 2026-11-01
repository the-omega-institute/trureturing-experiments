# Actual polygon ear slide

This companion checks the diagonal-set correspondence in B-37.2, equation (37.8), of the [FIB polygon source](https://github.com/the-omega-institute/trureturing/blob/b9087e51948fa40e487ce8039c2b126d0dff7a66/docs/develop/theory/AURIC_FIB_ATOM_INTRINSIC_CLASSIFICATION_AND_POLYGON_GEOMETRY.md). The source describes the full conditioned fibres: forced ears P, forbidden ears Q, and an exact observed ear trace, including absent observed ears. The local deletion correspondence is also described in Leader–Tan, *Some results on intersecting families of triangulations*, [arXiv:2610.39116v1](https://arxiv.org/abs/2610.39116), Theorem 2's proof. Neither reference is used as executable code.

From this private repository root, with Python 3.12 or newer and the standard library:

```sh
python3 experiments/fib-polygon-ear-slide-24755/check_actual_slide.py --output /tmp/fib-polygon-ear-slide-results.json
cmp experiments/fib-polygon-ear-slide-24755/results.json /tmp/fib-polygon-ear-slide-results.json
```

Assertions must remain enabled. A nonzero exit fails reproduction. The output is deterministic and uses exact integer vertex labels and finite sets.

For a cyclic block starting at a, let s be the ear centred at a+1 and t the ear centred at a+2. The map removes t, swaps these two vertices in every remaining chord, and inserts s. Its inverse removes s, performs the same vertex swap, and inserts t. Source and target triangulations are generated separately by root-triangle recursion; the target is filtered from that independent enumeration.

The default run checks:

- Every triangulation and every cyclic block placement for n=5 through 11: 21,682 source applications, with internality, noncrossing, maximality, equality with the entire ear-containing target, both inverse laws, and remote-ear membership.
- Every pair of distinct chords avoiding deleted vertex 2 for n=5 through 30: 509,730 crossing comparisons, plus the diagonal side exception and exclusion of chords incident to vertex 2.
- Eight literal conditioned placements at n=14, 15 and 16, covering (|P|,|Q|)=(2,2),(3,3),(2,4),(4,2), with the fourth local ear observed and unobserved in each case: 35,622 source completions and their independently enumerated target completions. Every remote exact-trace group, including its negative observations, is compared as an actual set of triangulations. Prescribed compatible ears are deleted for enumeration and restored as their original chords; forbidden ears are filtered on each side.
- The pentagon pairs {13,03} to {02,03} and {13,14} to {02,24}.

`results.json` contains the ranges, individual placement data and counts. These finite checks do not prove the theorem for all n, all P/Q/J/R, or arbitrary conditioned placements. In particular, the conditioned tests cover only the eight listed placements; the all-rotation enumeration has the smaller n=5 through 11 range. No Lean compilation, kernel proof, formal coverage or physical FIB realization follows from these experiments.

The program is repository-authored and imports no third-party code. It is covered by this repository's Apache-2.0 license.

## Literal forced-ear fans

`check_forced_fan.py` checks a second construction on the original cyclic vertex labels. For an independent set P of ear centres and a root h outside P, F(P,h) consists of the ear diagonals centred in P together with every internal diagonal from h to a vertex outside P. Independence means that no two centres are cyclic neighbours. The construction does not introduce or relabel a smaller polygon.

Reproduce its deterministic result with assertions enabled:

```sh
python3 experiments/fib-polygon-ear-slide-24755/check_forced_fan.py --output /tmp/fib-forced-fan-results.json
cmp experiments/fib-polygon-ear-slide-24755/forced_fan_results.json /tmp/fib-forced-fan-results.json
```

The default range is n=14 through 18. The program enumerates every independent P of sizes 2, 3 and 4 and every root outside P: 5,976 forced sets and 78,042 literal fans. It verifies internality, the n-3 diagonal count, noncrossing, and maximality against every original internal diagonal, using endpoint interleaving to compute crossings. The actual ear centres lie in P together with the two neighbours of h. For each of the 40,729 roots whose neighbours avoid P, exactly n-|P|-3 fan diagonals cross the ear centred at h; the minimum over the checked range is 7.

For every independent union K and ordered partition into disjoint P,Q with sizes (2,2), (3,3), (2,4) or (4,2), the program also checks two consequences. First, a root outside P whose neighbours avoid Q exists, and its fan forces P and avoids Q: 153,006 conditioned placements. Second, for every h outside K and its cyclic neighbours, the fan avoids Q and has no ear centre in E_K other than possibly the neighbours of h, where E_K excludes K and its cyclic neighbours: 413,250 strict witnesses. It is therefore empty on every observation H contained in E_K that contains h and excludes its neighbours. More than one diagonal crosses the ear at h, excluding the necessary unique-crossing property of a single-ear flip output. The program verifies that crossing count; it does not implement or prove existence of a flip.

`forced_fan_results.json` records the per-order counts and checked ranges. These are finite checks of the literal construction, not a theorem for arbitrary n, a Lean proof, or formal coverage of the source inequality. The program is repository-authored, imports no third-party code, and uses only the standard library under this repository's Apache-2.0 license.
