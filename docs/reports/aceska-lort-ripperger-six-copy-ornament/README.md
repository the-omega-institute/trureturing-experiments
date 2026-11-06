# Six-copy ornaments: exact construction checks

Source: Roza Aceska, McKailyn Lort and Allison Ripperger, *Ornaments and Difference Distance Magic Oriented Graphs*, arXiv:2606.16001v3, Section 1.2 and Section 4. This report uses v3's broad attachment definition: arbitrary oriented cross edges, preserving the old graph and its specified DDM labeling.

Run `python3 check.py`. The checked output is `check-output.txt`, ending in `ALL_CHECKS_OK`. Python 3's standard library suffices. Every computation uses integers, with no floating point or random inputs.

The fidelity checks reproduce the zero weight vector of the labeled wheel in Figure 1(a), transcribing its eight TikZ edges, and Remark 1's zero total cycle weight and nonzero individual vertex weights. The cycle checks cover all permutations for orders 3 through 7 at offsets 0 and 13: 11,820 labelings.

The construction checks cover every labeled oriented simple graph on one through four vertices: 1, 3, 27 and 729 cases, respectively (760 total). They also include an inward four-leaf star and a five-vertex transitive tournament. For each graph, three old orders are checked: the computed threshold, threshold plus 1, and threshold plus 17. The old graph is the nontrivial source wheel on labels 1 through 5 with isolated vertices elsewhere.

The checks verify complementary labels on six identically oriented copies, the two integer moment conditions, explicit ternary-generator decompositions of both affine coefficients, distinct old labels, actual cross-edge orientations, zero vertex weights and a cross edge at every added vertex. Sparse rows represent the cross-edge matrix; every unselected old vertex has a zero row by construction. This permits exact verification even though the elementary zero-contribution gadgets give large old-order thresholds. The largest threshold among the exhaustive small graphs is 28,856,472,338,804,772,962,515; no minimality or efficient threshold is claimed.

These finite computations check an explicit general construction and source fidelity. They are not a proof of the unbounded existence assertion or a Lean kernel check. The proposed settlement's main content is a general integer-lattice and label-realization argument, not positive finite enumeration. The original hypothesis's labeling is not imposed on the new copies, whose labels are chosen afresh.

The checker is original code; no source-paper program is copied. Mathematical formulas and the wheel's edge data are attributed above. The source is available at https://arxiv.org/abs/2606.16001v3.
