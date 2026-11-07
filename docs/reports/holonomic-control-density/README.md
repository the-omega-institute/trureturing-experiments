# Exact algebra checks for holonomic-control density

Target: Paolo Zanardi, *Universal holonomic control of algebras of observables*, arXiv:2609.28662v2, Appendix F, Conjecture 1 (Generic attainment).

Run with Python 3 and SymPy:

```sh
python3 docs/reports/holonomic-control-density/check.py
```

The program uses integer, rational and symbolic arithmetic throughout. The stored output is `output.txt`; successful execution ends with `ALL_CHECKS_OK`.

The fidelity check constructs an actual rational unitary two-band family by a Cayley transform and reproduces the paper's curvature formula and zero total curvature trace. The candidate check verifies the projected-commutator span for all 57 ordered positive band-size tuples of total dimension 2 through 6 having at least two bands. It also checks the elementary generator identities and the quadratic curvature scaling at a zero of the local field. A rational block-unitary transport model checks that a finite curvature determinant retains a nonzero leading coefficient and is nonzero at a specified rational parameter.

These finite checks support the algebra used in the proposed proof. The transport model is explicitly prescribed and is not claimed to be the Berry–Wilczek–Zee transport of the constructed local family. The program does not prove a statement about all band sizes, analytic dependence of actual parallel transport, Ambrose–Singer, or density in a smooth mapping space. Those remain mathematical proof obligations.

The program is independently written for this candidate. It does not execute or copy the source paper's numerical ancillary code.
