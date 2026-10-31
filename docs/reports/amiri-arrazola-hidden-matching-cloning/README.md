# Hidden-matching cloning: exact preregistration checks

Source: Ryan Amiri and Juan Miguel Arrazola, *Quantum money with nearly optimal
error tolerance*, arXiv:1610.06345v1, equations (9), (20), and (21).
DOI: 10.1103/PhysRevA.95.062334. The source conjectures the average marginal
fidelity bound 1/2 + 1/n for positive even dimensions.

Run `python3 check.py` with Python 3 and SymPy. `check.out` contains its output.
The program uses integer, rational, and exact symbolic arithmetic throughout.
No numerical eigensolver or floating point tolerance supplies evidence.

The checks reproduce the source's fidelity/error identity at n=4,6, using all
bit strings and a maximal disjoint family of perfect matchings. The average
marginal-fidelity bounds for n=4,6,8,10,12,14 are respectively 3/4, 2/3, 5/8,
3/5, 7/12, 4/7. All dimensions 2 through 14 are checked, including odd dimensions
as a test of the proposed stronger argument.

For n=2 through 8 the fourth sign moment is enumerated independently over every
source bit string. The resulting operator is decomposed into every singleton
and every three-distinct parity block, whose characteristic polynomials are
checked exactly. For n=9 through 14 the moment uses the even-multiplicity rule.
The 3 by 3 zero-sum block and 4 by 4 uniform block characteristic polynomials
are also checked symbolically in n. The proposed top eigenvectors, norms,
orthogonality, and identity input partial trace are checked at every tested
n. At n=2 through 6 the two output fidelities of their Kraus channel are
computed for every bit string and equal 1/2 + 1/n individually.

These are checks of an ordinary proof outline, not a proof for all dimensions
or a Lean settlement. The upper bound's proposed uniform argument uses the
fourth sign moment and invariant parity sectors. The branch contains no Lean
proof of the assertion or its negation.
