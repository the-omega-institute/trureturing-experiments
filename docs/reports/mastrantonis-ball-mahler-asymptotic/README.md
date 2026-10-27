# Exact checks for the Euclidean-ball Mahler asymptotic

Source: Vlassis Mastrantonis, *L^p-Legendre Transforms and Mahler integrals: Asymptotics and the Fokker–Planck heat flow*, arXiv:2411.10439v2, Conjectures 5.1–5.2 and Corollary 3.8.

Run with Python 3 and SymPy:

```sh
python3 docs/reports/mastrantonis-ball-mahler-asymptotic/check.py
```

The program reproduces the source's exact half-integer Bessel-K integral (TeX label KintAsympt1) for every odd n from 1 through 19 and m from 0 through 10. It independently evaluates the Gaussian integral underlying the source's convex-quadratic Mahler formula for p in {1/4, 1/2, 1, 2, 4} and n from 1 through 8.

It checks ball marginal coefficients through degree 24 in every dimension 1 through 12 against the source's normalized Bessel-I series, including direct polynomial integration in odd dimensions. It also checks Corollary 3.8's radial normalization on the same dimensions and five exact parameters. Finally it verifies, symbolically for p>0, the critical point, curvatures, exponential base and prefactor of the proposed two-stage Laplace argument.

Every comparison uses exact integers, rationals or symbolic expressions; no floating-point approximation supplies evidence. The stored output is check.out. This is an algebraic check of the expected settlement and source fidelity. It does not prove the dimension-asymptotic assertion: uniform inner Laplace estimates and outer tail bounds require analytic proofs. The proposed stronger equivalent is not an exact identity at finite n. No proof in Lean is included.
