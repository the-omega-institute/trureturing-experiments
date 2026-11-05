# Carnero Bravo Question 22: exact finite diagnostics

Source: Andrés Carnero Bravo, *On independence complexes of graph products*,
arXiv:2505.06457v2, Question 22 and Theorem 6.
<https://arxiv.org/abs/2505.06457v2>

Run with Python 3 (standard library only):

```sh
python3 check.py
```

The stored output is `check.out`. The program uses integer masks and integer
row elimination computing ranks over the rationals; it uses no floating point.

- Source fidelity: two individually verified elementary simplicial collapses
  take the independence complex of the five-vertex path to a four-cycle,
  reproducing Theorem 6's homotopy type S¹ for this instance.
- All 3,135 choices of labelled outer graphs on 1–4 vertices, labelled inner
  graphs on 1–3 vertices, and substituted vertices satisfy the deletion,
  cone-intersection, substitution-union and graph-link identities. All 159,322
  induced restrictions satisfy the corresponding smaller-substitution identity.
- Rational Betti additivity agrees with the proposed component splitting in
  every case with a nonempty link realization. Empty links are checked separately
  by the same exact face identities (451 cases).
- For 2 ≤ m,n ≤ 4, the independence complex of K_m ∘ (n K₂) is exactly a disjoint
  union of m boundaries of n-dimensional cross-polytopes. Its ordinary Betti
  vector has b₀ = b_(n−1) = m and all other entries zero. These are examples
  within the expected affirmative conclusion, and obstruct strengthening W to SW.

These are finite combinatorial, collapse and rational-homology diagnostics.
They do not certify the universal hereditary wedge-closure assertion. Rational
homology alone does not determine a homotopy type. The universal topological
argument and its Lean formalization remain separate obligations.

Original check program: MIT. Original output and README: CC0-1.0.
