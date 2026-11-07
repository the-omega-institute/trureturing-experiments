# Bipartite grid codes and rectangle transvections

Source: Victor V. Albert, *Beyond transversality: structure of Clifford circuits for CSS codes*, arXiv:2608.05688v2, §VI.2, equation (34) and its following conjecture. All arithmetic in `check.py` is exact binary arithmetic using Python integers.

Run `python3 docs/reports/albert2608full-grid-code-clifford/check.py` from the repository root. The expected final line is `ALL_CHECKS_OK`; the complete output is in `output.txt`.

The check reproduces the fourteen unordered admissible grids with `ab <= 120`, their code ranks, logical dimensions and distance four. It computes the entire transvection closure for the smallest source instance `(4,6)`: all 65,535 nonzero vectors of its 16-dimensional logical phase space occur. The corresponding symplectic group order is 59,980,383,884,075,203,672,726,385,914,533,642,240,000. The transposed `(6,4)` instance checks the coordinate symmetry. The `(4,4)` instance is an additional control outside the source's stated parity range.

For all fourteen grids, the program checks every rectangle's identity `I + R + Q + RQ = w w^T` on every physical coordinate vector; verifies that row and column transpositions preserve the classical code; checks the nondegenerate logical Gram matrix, spanning, connected pairing graph and an independent orthogonal sum triple. It also tests the saturation enlargement construction in all 240 configurations of a four-dimensional saturated summand and one additional hyperbolic plane.

These finite checks do not prove the infinite-family conjecture. Full transvection closure is enumerated only for the smallest instance and the control. The remaining thirteen source instances receive structural checks, not a fresh computation of their whole logical groups. The proposed all-parameter argument requires proofs of the quotient identification, physical-layer construction and abstract saturation lemma. No Lean proof of the conjecture or its negation is included.

The program is independently authored from the source's mathematical definitions; no third-party program source is copied or imported.
