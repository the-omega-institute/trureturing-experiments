# Finite building code duality: exact checks

`check.py` implements the chamber-indicator definitions in Nolan J. Coble,
*New Directions From Old Codes* (2025), Appendix B, Definitions B.7 and B.9,
Question B.11 and Proposition B.12. Primary source:
https://nolanjcoble.com/assets/files/Coble_Dissertation.pdf
DOI: 10.13016/vsd7-e6no.

Run with Python 3.12 or newer, without dependencies:

    python3 docs/reports/coble-2025-question-b11-finite-building-code-duality/check.py

`output.txt` is the resulting stdout. Vectors are integer bit masks, binary
row operations use XOR, and projective incidence and normalization use integer
arithmetic modulo 2 or 3. No floating point is used.

The test set consists of 13 explicitly constructed buildings: rank zero;
rank-one buildings with 2, 3, 4 and 7 chambers; products with panel sizes
(2,2), (3,3), (2,3,4), (2,3,5); thin types A2 and A3; and the full-flag
buildings of GL3(F2) and GL3(F3). It exhausts all subapartments and residues
of every rank in each listed building and checks all orders from -1 to m.
It reproduces the source's Fano and K3,3 examples, its orthogonality inclusion,
its endpoint identities and the thin-code specialization. Equality follows
in each instance from exact orthogonality and complementary dimensions.

It also constructs the proposed triangular apartment basis and checks its
simultaneous panel kernels, all residue kernels and all descent thresholds.
Products with unequal panel sizes check that the computation does not rely
on one common thickness. These finite results do not prove the assertion
for arbitrary finite buildings and do not certify the building axioms in Lean.
The constructions and the row-reduction implementation are written for this
report; no third-party program source is copied.
