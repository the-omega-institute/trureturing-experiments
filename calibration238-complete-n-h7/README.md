# Complete normalized seven-coordinate acquisition certificate

This directory supplies the exact dataset and independent arithmetic for the
normalized slice

$$
1\ge b_1\ge b_2=b_3\ge b_4\ge b_5\ge b_6\ge b_7\ge e\ge0,
\qquad 0\le b_3+b_4+b_5+b_6+b_7-8e\le1-b_1.
$$

The associated ordinary theory is *Self-calibrating acquisition resources XIX*,
Sections 341–347. The coefficient inverse is an ordinary general proof; the finite
certificate establishes the full normalized equality-slice classification. The
scientific repository is `the-omega-institute/trureturing-experiments`, and the
relative scientific directory is `calibration238-complete-n-h7/`.

`H7-certificate.json` is the entire Research161 H7 dataset, schema
`research161-H7-complete-normalized-common-allocation-v1`, with 2,579,533 bytes and
SHA256 `621bed2ac570875c2bcb96796338e753ba045c544aa30fb8e1af8702cc94581d`.
Its original certificate bytes are retained exactly. The data comprise ten domain
inequalities, 471 full preselected-class supporting scores, all 1,606 vertices of
their intersection, and one complete normalized stopped allocation at each vertex.
No individual cut or vertex is attributed as novel here.

`example-15-over-64.json` gives the explicit fourteen-schedule allocation and the
complete supporting score for prefix `(43,34,34,29,27,25,22)/64`. It is a member of
the classification. Every positive ordered vanishing tail with `a8 <= 15/64` is
allowed in the original marked-target interpretation. The allocation and supporting
score are displayed in Section 346 of the associated theory; their arithmetic is
independently reconstructed by the verifier.

Run with Python 3.9 or later, using the standard library only, from this directory:

```sh
python3 verify_h7.py --certificate H7-certificate.json --example example-15-over-64.json --output results.json
```

Every file path is an explicit argument. Verification is offline and uses no
optimizer, certificate generator, personal path, project import or hidden data.
An invalid identity or inequality raises an exception and exits nonzero; a complete
successful run writes the mathematical counts to the explicit output path and
prints that same result. `results.json` is the result of the complete calculation,
not a sampled test or an optimization transcript.

The precise exact scope is as follows:

- Reconstruct the 382 event rows on all 128 words and eight labels. For each score,
  sum its nonnegative weights into the whole reward array. Compute the 21 pair
  supports by taking all eight label maxima separately on each of the four
  transcript cells. Check the coefficient quotient, constant demand, total weight,
  and positive primitive normalization of each of the 471 cuts. Singletons,
  constants, repetitions and reversal are covered by information refinement of a
  pair, as proved in the theory.
- Solve all 120 seven-row domain bases over the rationals, obtaining the twelve
  stated initial vertices. Insert every cut in order. An edge is recognized by
  rank six of the intersection of full active sets. Retain inside/on vertices and
  add every strict crossing-edge intersection. Coordinates use primitive integer
  homogeneous vectors; ranks use integer elimination. There are 21,337 crossings,
  a peak of 1,752 vertices, and 1,606 final vertices with 8,363 final edges. Compare
  the complete generated vertex set and active masks to the certificate, test all
  481 inequalities at every final vertex, and prove final affine rank seven.
  Exactly 251 vertices have exactly seven active rows. Full final dimension also
  ensures that every intermediate containing polytope was full dimensional. The
  face and clipping proof in the theory makes this a complete enumeration rather
  than mere membership checks of supplied points.
- Bind every lift to its vertex. Check all three counterfactual cells separately
  against one schedule mass, all labels' nonnegativity, and the global mass-one
  budget. Evaluate every word and all independent base, correct, union and positive
  inequalities literally. Totals are 22,027 schedule rows, 66,081 cell vectors,
  205,568 word laws and 613,492 independent event inequalities.
- Reconstruct the displayed example's 22-row score and all 21 four-cell supports,
  its 14 schedules and 42 cell vectors, all 128 laws and 382 event rows, exact
  correct counts, outside maxima, base deficits, positive law and stopped
  coefficient identity `c=H`. The lower cut is
  `b1+4*b3+5*(b4+b5+b6+b7) <= 1+42*e`, with exact optimum `15/64`.
- Enumerate all 820 deterministic schedule/table descriptions on three bits and
  four labels for the scalar-budget obstruction. Exactly 354 avoid label zero on
  every forbidden word, and their maximum two-word zero score is one, whereas the
  prescribed law requires `3/2`. This finite illustration supplements the ordinary
  all-schedule proof; it does not establish an original profile separator.

The closed real polytope follows from its complete vertex hull and convex mixtures
of these entire allocations. A single seed fixes a branch, schedule and every
counterfactual quantile before replies. Original-source conclusions use the
separate class-feasibility transport and finite-pattern theorem, not just the
finite arithmetic. There is no unequal-`b2,b3`, all-mass, other-horizon, facet
minimality, universal N/O, M1, nonatomic-target or current-kernel conclusion.

The literature used for the ordinary proof has the following exact roles.
Blackwell, *Comparison of Experiments* (1951), Sections 2–3 and Theorem 2, treats
whole attainable finite-dimensional loss vectors for finite families of Borel
measures and bounded action/loss sets. Its converse in that paper is binary.
Kaniewski–Lee–de Wolf, arXiv:1411.7280v1, Definitions 6/8 and Theorem 9, treats
nonnegative scalar expectations at maximum supported tree depth; its coefficient
normalization rescales scalar output and supplies no equal vector-cell budgets.
Kaibel, arXiv:1001.2161v1, Theorems 8/13/15 and Section 1.5, supplies the mature
halfspace, face, vertex and bounded-hull facts. These do not supply this dataset.
The full primary texts were accessible at their cited versions; relevant
hypotheses were inspected separately from this offline calculation.
