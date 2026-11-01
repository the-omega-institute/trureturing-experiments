# Exact checks for the orbital majorization-least conjecture

Source: Aliverti-Piuri, Liebert and Schilling, *Maximal Orbital Entanglement from Non-Interacting Fermions*, [arXiv:2610.32974v1](https://arxiv.org/abs/2610.32974v1), Conjecture 1. The program uses Python 3 standard-library integers and `fractions.Fraction`; it uses no floating point or numerical eigensolver.

Run `python3 docs/reports/apls-conjecture1-majorization-least/check.py` from the experiments repository. The deterministic output is stored in `output.txt`.

The fidelity checks reproduce the source's dimension-two half-occupation optimum and its dimension-three `(2,1)` formula. Exact idempotence, trace and determinant certify the first matrix's spectrum; exact sum/difference identities certify the second matrix's spectrum.

The main check covers every admissible partition for dimensions two and three and every decreasing occupation vector on the grid `{0,1/8,...,1}`: 210 occupation inputs and 375 continuous feasible polytopes. The feasible polytopes are the complete Schur-Horn diagonal polytopes and the Cauchy interlacing rectangle stated in the source. Coordinate-equality hyperplanes divide each polytope into chambers where sorted partial sums are affine. Exact intersections of boundary/chamber hyperplanes include every chamber vertex, including lower-dimensional degenerate cases. Verifying the majorization inequalities at those vertices verifies every feasible point for each tested occupation input.

The program also checks membership in the finite set of repeated subset averages and exact supporting-gradient inequalities for the objectives `-sum(x_j^2)` and `-sum(x_j^2+x_j^4)`. Their gradients at the proposed optimum support the entire polytope. Linearity in the objective parameter gives the same certificate for every member of their convex homotopy.

Result: 2,587 chamber vertices, 7,644 partial-sum comparisons and 5,174 supporting-gradient checks passed, including endpoint and repeated occupations. This finite evidence does not prove the unrestricted conjecture, test dimensions above three, establish spectral convexity in higher dimensions, or validate the spectral trace derivative. No binary entropy is approximated numerically. The planned positive settlement is a general analytic argument; these checks are not its main content.
