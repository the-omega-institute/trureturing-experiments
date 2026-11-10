# Deterministic-arm necessity: exact preregistration checks

Source: Ivan Viakhirev, Kirill Borodin and Grach Mkrtchian,
*Sharp Bounds on the Winner's Bias under Independent and Nonnegatively Correlated Noise*,
arXiv:2610.31948v1, Theorem 3, Section 7 and Appendix D.5.

Run `python3 check.py` with Python 3.12 or newer and SymPy 1.14 or newer.
The program uses rational arithmetic and symbolic polynomial identities only.
Its stored output is `output.txt`; success ends with `ALL_CHECKS_OK`.

The fidelity checks reproduce the source's cascade expectation 1/2, total
variance M/(2(M-1)), and sharp squared ratio (M-1)/(2M) for M=2 through 10.
The Appendix E sum-of-squares identity is checked symbolically for arbitrary
parameters; the survival telescoping identity is checked symbolically through
five arms.

For M=2 through 5 the program enumerates unordered tuples from nine
nondegenerate centered two-point laws (positive probability 1/3, 1/2 or 2/3;
positive value 1/2, 1 or 2), together with the zero law. It checks 2,992 tuples,
including strictness for all 1,992 tuples with no zero arm, and the first-positive
partition and Cauchy–Schwarz coefficient bounds. Positive-variance equalities
found on this grid all contain a zero arm. Scaled, permuted source cascades
separately check actual equalities for every M=2 through 5.

Four finite centered laws test the explicit positive/negative product mixture.
The marginal reconstruction, total mass, average conditional variance, expected
maximum and positive mass of all-nondegenerate product components are exact.

These finite checks support the proposed general argument. They neither
classify all probability laws nor prove the general conjecture. The proposed
settlement needs an argument transferring strictness through measurable
centered product mixtures, beyond the source's supremum-preserving reduction.
