# Zimmermann–Stoye codimension-one Euclidean Stiefel curvature

Source: R. Zimmermann and J. Stoye, *High curvature means low-rank: On the sectional curvature of Grassmann and Stiefel manifolds and the underlying matrix trace inequalities*, [arXiv:2403.01879v3](https://arxiv.org/abs/2403.01879v3), Conjecture 1, equation (26), and Theorem 10(5). The metric is the ambient Frobenius metric, with no factor 1/2.

Run with Python 3 and SymPy:

```sh
python3 docs/reports/zimmermann-stoye-stiefel-curvature/check.py
```

`output.txt` contains the output. Every calculation uses exact integers, rationals, or symbolic algebraic numbers; no floating point arithmetic is used.

The fidelity checks reproduce the displayed curvature 1/2 of Theorem 10(5) on St(3,2) and the curvature -1/2 of the endpoint blocks in equation (28). For unrestricted skew blocks in dimensions p=2,3,4,5, the program checks equation (26) against the embedded Gauss numerator, the normalized-row expansion, the cancellation in the coupling term, and its Cauchy sum-of-squares identity and half-Frobenius budgets. It also checks the scalar two-square upper-deficit identity, exact endpoint blocks for p=2,...,8, and every coordinate tangent 2-plane for p=2,...,5.

These are diagnostics for the proposed dimension-independent argument. The finitely many symbolic matrix dimensions and coordinate planes do not prove the universal conjecture. Ge's dimension-independent skew-commutator bound and the orthogonal normalization are mathematical ingredients of the proposed proof, not certified by these finite checks. No Lean proof is included.
