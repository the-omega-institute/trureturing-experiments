# Balanced-vector recognition under binary encoding

The source is Johnston, Plosker and Varona, *Enumeration of Laplacian integral and {-1,0,1}-diagonalizable graphs*, arXiv:2607.06336v1, Definition 2, Section 8's fourth question, and Appendix A.1. The expected answer is NP-completeness for explicitly listed positive binary integers. No Lean proof is included.

Run `python3 check.py` with Python 3; no third-party packages are needed. All calculations use exact integers. `output.txt` is the resulting output, ending in `ALL_CHECKS_OK`.

The fidelity check reproduces the balanced-partition counts for total sums 1 through 9 in the source's table: 1, 2, 2, 4, 3, 7, 6, 13, 12. A direct enumeration of ternary annihilators independently checks the generating-function Gram-matrix recognition oracle on 125 vectors. That oracle computes the rank of the sum of outer products of every ternary annihilator, using exact coefficient extraction and polynomial factor deletion.

The construction checks cover the appended-coordinate criterion, 234 balanced completions with full-rank ternary certificates and all auxiliary entries above the target, 72 complete bounded-two subset-sum reductions with independently recognized final vectors (also after sorting), 282 no-carry reductions from ordinary subset sum, and three composite reductions containing both outcomes. A binary-size example uses 127-bit inputs and yields 390 completion coordinates of at most 129 bits.

These finite results check the formulas and expose small-case failures. They do not prove the universal reduction, NP-completeness, the single-tape running-time bounds, or novelty. The main proposed argument is the general balanced-completion and saturation construction; the finite enumeration is verification evidence only. The pseudopolynomial recognition oracle also distinguishes the binary and unary encodings: an algorithm polynomial in the total coordinate sum need not be polynomial in its binary length.
