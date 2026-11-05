# Jeong's concave Jordan orbit inequality

This report checks finite instances of Juyoung Jeong, *The Schatten p-quasinorm on Euclidean Jordan algebras*, AIMS Mathematics 9(2), 5028–5037, DOI [10.3934/math.2024244](https://doi.org/10.3934/math.2024244), published Conjecture 4.2. It contains no Lean proof and does not establish the universal conjecture.

Run with Python 3, using only its standard library:

```sh
python3 docs/reports/jordan-concave-orbit-inequality/check.py
```

`check.out` is the complete successful output. All evidence uses integers and `fractions.Fraction`; no floating point or tolerance is used.

- Fidelity: reproduce the equality in Theorem 3.5 for two orthogonal primitive idempotents at p = 1/2, obtaining 4 = 2(1 + 1).
- Check the orbit inequality using identity automorphisms for every pair of vectors with entries in {−2,−1,0,1,2}, in componentwise Jordan algebras of ranks 1–3, for eight specified increasing nonnegative concave functions: 130,200 checks.
- Check every ordered pair from the 53 real symmetric 2×2 matrices with entries in {−2,−1,0,1,2} and rational eigenvalues, retaining the 1,561 pairs whose sum also has rational eigenvalues. There are 224 noncommuting pairs and 437 singular sums. The eight functions give 12,488 orbit certificates. The 1,248 pairs whose sum has irrational eigenvalues are explicitly excluded.
- Check the scalar contraction Jensen inequality at 81 rational parameter pairs for each function: 648 checks.

For the matrix checks, let α and β be the decreasing eigenvalue pairs of f(|A|) and f(|B|), and let c be those of f(|A+B|). The checker chooses z in

\[
[\max(\alpha_1+\beta_2,\alpha_2+\beta_1),\alpha_1+\beta_1]
\cap[c_1,\alpha_1+\alpha_2+\beta_1+\beta_2-c_2].
\]

It constructs the squared cosine t of a real rotation and verifies 0 ≤ t ≤ 1, the trace, and the determinant of diag(α) + R diag(β) Rᵀ exactly. The resulting eigenvalues z and total − z dominate c componentwise. Real orthogonal diagonalization aligns this sum with f(|A+B|), providing the required orbit witnesses. Squared off-diagonal entries avoid irrational arithmetic. The trace inequality from Theorem 3.4 is also checked independently on these instances.

Functions with positive value at zero, the function 1 on (0,∞) and 0 at zero, and its translate by 1 are included. The checks therefore retain the singular and endpoint-discontinuous cases. Exceptional algebras, irrational-spectrum matrix pairs and arbitrary functions are outside these finite checks; the proposed universal settlement requires a general argument.

The program and report are original experiment material. The source article is cited rather than redistributed.
