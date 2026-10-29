# Pure-state coevent separation check

This report gives an exact finite check for Petros Wallden, “Distinguishing Initial State-Vectors from Each Other in Histories Formulations and the PBR Argument”, arXiv:1211.2084v2, Appendix. The Python 3 program uses `Fraction` pairs for exact real and imaginary coefficients (with a common `1/sqrt(2)` factor) and SymPy exact identities; it uses no floating point.

It reproduces the source protocol’s minimal null events and primitive supports for three normalized rational Pythagorean overlaps, checks the phase-regularized protocol’s same disjoint classifications, and verifies that adjoining one or two zero branch atoms cannot occur in a primitive support. The finite checks support the proposed analytic route; they do not by themselves prove the universal claim.
