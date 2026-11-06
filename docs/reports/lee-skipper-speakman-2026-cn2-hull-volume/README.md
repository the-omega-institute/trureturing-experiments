# Monomial graph hull with two nonzero lower bounds

`check.py` checks the proposed finite inequality system and volume formula answering the Section 4 Outlook question of Jon Lee, Daphne Skipper and Emily Speakman, [arXiv:2605.01493v1](https://arxiv.org/abs/2605.01493v1). It also reproduces the source's Section 2 Theorem 1 inequalities and Section 3 Theorem 2 volume formula for the one-nonzero-lower-bound special case.

Run with Python 3, SymPy and SciPy:

```sh
python3 check.py
```

`check-output.txt` is the program's output. All evidence for facet coverage and volume uses integers or exact rational arithmetic. For each tested box, the program exhausts the affine hyperplanes through every dimension-sized subset of lifted corners, retains exactly the supporting facets, and checks that every such facet occurs in the nominated system and every nominated inequality is valid at every corner. It computes volume independently by an exact pulling triangulation and rational determinants. Thus each tested inequality system describes exactly that corner hull.

The 40 boxes comprise all nine pairs of normalized lower bounds in {0, 1/3, 2/3}, for each n in {2, 3, 4}; three seeded random rational boxes in each dimension; the n=4 half-bound example with volume 31/640; and one additional one-nonzero-lower-bound box in each dimension. Twelve boxes reproduce the source's complete C1 hull and volume results. Symbolic algebra also checks the envelope-integral expression against the proposed volume formula for every integer n from 2 through 20.

As a separate numerical diagnostic, SciPy ConvexHull computes each corner hull. Every nondegenerate simplicial facet is reconstructed exactly from its rational vertices, and the recovered set must equal the exhaustive exact facet set. Rank-deficient boundary simplices emitted by Qhull on flat facets are discarded. Its floating-point volume is compared to the exact volume with relative/absolute tolerance 1e-8. Floating-point agreement is not used as proof evidence.

The source's graph-to-corner-hull equivalence is an input to these finite tests. The tests are neither a general-dimensional proof nor a Lean proof. The proposed all-dimensional settlement requires the finite common-distribution construction and elementary envelope integration described in the preregistration.
