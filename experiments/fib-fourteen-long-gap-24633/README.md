# Fourteen-gon long-gap candidate

This companion reproduces the exact finite computation for the proposed long-gap matching-core bound in §39 of `docs/develop/theory/AURIC_FIB_ATOM_INTRINSIC_CLASSIFICATION_AND_POLYGON_GEOMETRY.md` in the main trureturing repository. The mathematical domain, ear-deletion and inclusion-exclusion argument, arbitrary-family bridge, maximal-clique recurrence and full table belong to that manuscript. This directory supplies the necessary executable program and finite results, without a duplicate manuscript snapshot.

From this private repository root, with Python 3.12 or newer and no third-party dependencies:

```sh
python3 experiments/fib-fourteen-long-gap-24633/verify.py --results /tmp/fib-fourteen-results.json
cmp experiments/fib-fourteen-long-gap-24633/results.json /tmp/fib-fourteen-results.json
```

`verify.py` uses exact integers. It checks 843 independent cycle masks, 828 positive full-ear weights, and all 208012 actual triangulations counted with multiplicity; the direct root-triangle recursion must agree with every inclusion-exclusion weight. For each of all 15 pair partitions, the incremental maximal-clique recurrence and a separate pivoted enumeration must return identical entire maximal-clique sets. All three exact-trace core totals are checked against 4344. `results.json` records all fifteen rows, including vertex counts, all/spanning maximal-clique counts and maximum spanning weights. No floating-point arithmetic, solver, symmetry quotient or selected-completion sample is used. Assertions must remain enabled; a nonzero exit is a failed reproduction.

The supplied source is `mathematical-result.md`, artifact `6ae585c69b5dc1f88b8eccdc`, SHA-256 `1ffc20ecba6f58312dc7c11a750145ac6a172e6d6dd925c676dd7d6ea1c3aaa0`. The verifier is its complete Python block, with result serialization added; the enumeration, expected triples and assertions are preserved. Project parameter sources are B-32/B-33/B-37/B-38 in main commit `3143b5b6d2dcd2c00a95c3fb99ba5db59f498f1d` for B-32/B-33/B-37, and the manuscript's unchanged B-38 for its separate gap pattern. The finite program is repository-authored source supplied for preservation; it imports no third-party code.

These results are author reproduction of a supplied candidate, not independent mathematical acceptance. The candidate asserts A <= 4290 for all three active exact matching cores in any labelled strictly convex fourteen-gon with selected-centre gaps (3,1,1,1,1,1). It retains full ear sets and every diagonal completion; it neither requires family saturation nor asserts attainment of 4290. The A > 4344 residual would be empty, not an inhabited extremum. No n >= 15, Lean/kernel, formal admission, global novelty, physical realization or model diversity claim follows. Source production and exposition share GPT-family/repository exposure; prior provider memory is unknown.

The quoted source comparisons are Leader–Tan, *Some results on intersecting families of triangulations*, arXiv:2610.39116v1 §§2–4, with T_(14-r) = C_(12-r), and Talbot, *Intersecting Families of Separated Sets*, arXiv:math/0211314v1, Theorem 4. They supply neither this exact full-ear-weight optimization nor independent verification of it; the manuscript states their distinct scopes.
