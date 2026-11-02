# Continuous Jordan–von Neumann reconstruction: exact checks

Source: Alexandra Paicu, Dorian Popa and Mircea Dan Rus, “An Integral Formula for the Inhomogeneous Jordan–von Neumann Equation”, arXiv:2606.24565v1, Example 1; version of record, Mediterranean Journal of Mathematics 23, 207, Example 4.1, DOI [10.1007/s00009-026-03206-z](https://doi.org/10.1007/s00009-026-03206-z). The journal article is available under CC BY 4.0. The program independently implements the displayed mathematical formulas; no source code from the article is included.

Run with Python 3 and SymPy (verified with SymPy 1.14.0):

```sh
python3 check.py
```

All arithmetic is exact, using SymPy rationals and polynomials or Python `fractions.Fraction`. The stored result is in `check.out`.

The program reproduces the source identity Q(x⁴/12) = x²y². It solves the coefficient constraints for every polynomial cocycle of total degree at most 8 (45 coefficients, rank 37, nullity 8), and verifies each reconstruction identity on a basis of the full solution space. All the operators checked are linear in the cocycle, so the basis covers this polynomial space.

The candidate construction uses the triangular kernel ρ(z) = max(1−|z|,0), F(x) = ∫₋₁¹ρ(z)g(x,z)dz, b = ∫₀¹g(0,z)dz−∫₋₁⁰g(0,z)dz, a(t) = g(t,1)+g(t,−1)−2g(t,0), u(x) = −bx/2 + (1/2)∫₀ˣ(x−t)a(t)dt, and f = u−F/2. A negative control checks that the integral without the (x−t) factor fails the required second-derivative identity.

Three continuous functions q(x) = |x−c|, for c = −1/3, 0, 1/3, supply unbounded, nonsmooth cocycles g = Qq. Exact piecewise Simpson integration is valid here because the integrands are polynomials of degree at most two on the explicitly split intervals. On the half-integer grid from −3 to 3, the program checks 6,591 cocycle triples and 507 reconstruction pairs, including negative and zero integral endpoints. It independently verifies that the reconstructed f differs from q by a homogeneous quadratic.

These symbolic and finite checks support formula fidelity and detect sign or normalization errors. They do not establish the universal existence statement for arbitrary continuous cocycles. No Lean proof is provided.
