# Werner–Holevo operator-convex trace checks

`check.py` uses Python 3 and SymPy, with integer, rational and exact Gaussian-rational arithmetic. It never uses floating point as evidence. Run:

```sh
python3 docs/reports/muller-werner-holevo-operator-convex/check.py
```

`check.out` is the stored output; a successful run exits 0 and ends with `ALL_CHECKS_OK`.

The source is Markus Müller, *Convex Trace Functions on Quantum Channels and the Additivity Conjecture*, arXiv:0809.4060v3, Theorem 10 and its following unnumbered conjecture, DOI 10.1103/PhysRevA.79.052332. The check reproduces the source's two-factor pure-product and maximally entangled output spectra for dimensions 2 through 6.

For all 15 listed dimension tuples, the program compares every bilinear coefficient of the channel-purity identity against a separate partial-trace matrix-unit formula. The tuple set consists of one factor of dimension 2, 3 or 4, and every two- or three-factor tuple with local dimensions 2 or 3. Thus the identity check includes the full matrix algebra in each selected size, rather than a sample of pure states.

The state checks evaluate the purity bound and trace bounds for the square, shifted rational resolvents, the source's rational operator-convex generators and an endpoint-jump function, on exact pure and mixed examples. They include complex coefficients, unequal dimensions, three factors, maximally mixed states, product-state attainment and all distinct projectors of nonzero real ternary vectors on two qubits. Matrix inverses and ranks use exact domain arithmetic. A pure-product projection identity and rank check certify attainment in the selected cases.

The purity identity is credited to R. Alicki and M. Fannes, quant-ph/0407033v2, display labeled `ineq`. The proposed arbitrary-function argument additionally uses F. Hiai, *Matrix Analysis: Matrix Monotone Functions, Matrix Means, and Majorization*, DOI 10.4036/iis.2010.139, Theorem 2.5.2 and Corollary 2.5.4. Those analytic results are not verified by this script.

These finite checks test the formulas and the proposed route. They are not a proof of the conjecture for arbitrary tensor length, dimensions or operator-convex functions. No Lean proof is included. The program is authored for this report; no third-party source code is copied. SymPy is an external runtime dependency distributed under its own license.
