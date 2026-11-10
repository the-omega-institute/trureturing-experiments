# Normality under flag completion

Exact checks for Benedetti–Menara, *Diameter bounds for complexes that are close
to flag or close to normal*, arXiv:2611.08071v1, Question 9.1.

Run `python3 check.py`. Python 3.12 or newer suffices; there are no dependencies.
The program uses integer cardinalities, sets, exhaustive maximal-clique search,
and breadth/depth-first graph search. The stored output ends in `ALL_CHECKS_OK`.

The fidelity check reproduces the cross-polytope-minus-a-facet example in
Section 8.1 for dimensions d = 3, 4, 5, 6: normality, G₁, its flag completion,
its minimal nonfaces, dual diameter d, and 1-local diameter 2.

For each 6 ≤ N ≤ 24, the second check constructs a pure 3-complex on N + 4
vertices. Let a = N, z = N + 1, w = N + 2, b = N + 3. Its ordered facets are:

- zwij for each edge ij of the cycle 01, 12, …, (N−1)0;
- az01, az12, az23, az34;
- wb01, wb(N−1)0, wb(N−2)(N−1), …, wb34;
- ab01.

Every attachment meets the preceding complex in a nonempty union of triangles;
the program also independently checks every induced facet-star dual graph,
including the empty-face star. There are 2N + 3 facets. The exact f-vector is
(N + 4, 4N + 8, 5N + 8, 2N + 3).

The flag completion adds exactly the tetrahedron ab34 and triangles ab3, ab4.
It stays pure and its full dual graph stays connected. Its ab-star consists
of ab01 and ab34, with no dual edge, so it is not normal. Its f-vector is
(N + 4, 4N + 8, 5N + 10, 2N + 4). At N = 6 these vectors are respectively
(10, 32, 38, 15) and (10, 32, 40, 16).

For d = 4 and k = 2, G₂ forbids minimal nonface sizes 3 through 2, an empty
range expressly allowed in the source. No ghost vertices, empty graph, loss of
purity, dimension increase, or disconnected full dual graph causes the failure.

The range 6 ≤ N ≤ 24 is finite computational evidence. The argument for every
N ≥ 6 uses the same ridge attachments and the disconnected overlap of the a-
and b-neighborhood paths on the cycle; an unbounded family is not certified by
this finite run. These checks are not a Lean proof or a priority claim.
