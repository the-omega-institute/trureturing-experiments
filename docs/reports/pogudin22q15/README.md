# Pogudin Question 15: exact checks

Source: Gleb Pogudin, *Some Open Problems in Differential Algebra*, revision headed January 24, 2025, Section 4, Question 15. The source is the author's web PDF, not an arXiv paper:
https://www.lix.polytechnique.fr/~pogudin/files/da_problems.pdf

`source-pages-6-8.pdf` preserves the three source pages containing the integrating-factor discussion, the definitions of both constructions, and Question 15. Text extracted from this excerpt agrees with the corresponding original pages after whitespace normalization. The source PDF has SHA-256 `ef478da480ab3abaec2157b19e81623f747c6828c941900dc4f3e1d81ea8714e`.

Run with Python 3.12 or later and SymPy 1.14.0:

```sh
python3 docs/reports/pogudin22q15/check.py
```

The checked output is `check.out`. It ends in `ALL_CHECKS_OK`.

The program uses exact integer and rational symbolic algebra. It reproduces the source identities D(xx') = xx'' + (x')² and D((x')²/2 + x) = x'(x'' + 1). It also reproduces the first construction's six-dimensional image closure: the symbolic map lies in the stated six-term span and has Jacobian rank six at an exact integer point.

For every checked integer k from 2 through 12, it verifies

- p_k = xz + (k/2)(y² + 1), H_k = x^k(y² + 1), M_k = 2x^(k−1)y;
- D(H_k) = M_k p_k, total degree two, and the primitive linear-in-z criterion;
- exclusion from the first construction by the nonzero constant coefficient;
- exclusion from the second construction by the full coefficient matrix for every affine A(x,y,z) and every cubic C(x,y,z), not sampling their coefficients: its kernel forces A = 0.

The chosen k=2 polynomial also fails the source's equations for both components of the second Zariski closure. The bounded first-integral calculations find precisely span(1,H_k,H_k²) through total degree 2(k+2), for k=2 through 6. These bounded checks do not prove the full first-integral-ring description or an unbounded minimum-degree theorem.

The source leaves the admissible class of integrating functions unspecified and explicitly allows non-rational expressions. Its integrability discussion must not be described as restricted to polynomial multipliers and polynomial primitives. These experiments establish explicit polynomial witnesses and algebraic non-membership checks; they are not a Lean settlement or a completed source-fidelity judgment.
