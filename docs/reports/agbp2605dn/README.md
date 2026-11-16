# Exact checks for the inner-perimeter Dirichlet–Neumann question

Run `python3 check.py` (Python ≥3.12 with SymPy). The adjacent `check.out`
contains the output and ends with `ALL_CHECKS_OK`. Every comparison uses
integers, rational arithmetic or symbolic polynomial algebra; there is no
floating-point calculation or numerical eigensolver.

The fidelity check reproduces the ball case of the Steiner formula and
`W_(n−1)(B_r)=ω_n r` from arXiv:2605.25182v1, section 2 (TeX labels
`eq:steinerRN` and `eq:quer-ball`), for dimensions 2 through 6. Further checks
verify the radial integrals, capsule area and shell-volume constraints, and
rational radius brackets. For `L=k²`, the analytic estimate `2r/R³` lies in
`(31,32)`, `(1457,1458)` and `(1459,1460)` for `k=16,729,730`, respectively.
The matching annulus has an admissible trial quotient strictly between
20 and 21. The symbolic family estimate separates every `k≥729`.

For an additional spectral comparison, the annulus radial Neumann equation
is `b μ cos(μ(b−1))−sin(μ(b−1))=0`, with `b³=5/2`. Rational interval bounds
and alternating Taylor remainders give a root in `(387/100,388/100)`, hence a
first radial eigenvalue in `(149769/10000,150544/10000)`. Positivity of the
radial solution identifies it with the principal eigenvalue. The script
checks the endpoint signs, rather than certifying this PDE identification.

These calculations verify the algebra in the proposed analytic refutation.
They do not certify the Sobolev trace/density theorem, capsule coordinate
changes, the Cauchy–Schwarz argument for arbitrary Sobolev functions, smooth
approximation, or any Lean proof.
