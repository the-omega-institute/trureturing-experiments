# Singular Lorentz-Gram regularization

The target is the final regularization-independence clause of Question 9.1 in
Mihai Putinar and Prateek Kumar Vishwakarma, *Hyperbolic distance matrix completion*,
arXiv:2609.10403v1. This report supports a preregistration draft; it is not a Lean
proof and does not settle the separate question of an intrinsic singular selector.

Run with Python 3 and SymPy 1.14.0:

```sh
python3 check.py
```

The stored output is `check.out`. All evidence uses exact integers, rational
numbers, polynomial identities, matrix identities, or symbolic rational limits.
No floating point is used.

The fidelity test reproduces Theorem D's product-distance completion on a
three-vertex tree, including its Lorentz realization, determinant product and
inverse sparsity. The boundary test checks a symbolic family on the connected
chordal graph K5 minus edge 45. Its separator determinant is 9 epsilon^2/16;
the canonical missing entry is 25/16 - ab. Two rational regularizations have
the same prescribed singular limit and completed missing entries 25/16 and
34/25. Their difference is 81/400. Each limiting maximal clique has rank 3;
the five points are distinct and every prescribed off-diagonal Gram entry is
strictly greater than one.

The 245 rational instance checks supplement the symbolic identities. They
are not an enumeration-based proof about all regularizations. The proof
outline must additionally connect the inverse-sparsity or clique-transfer
characterization to the source's canonical completion in Lean. The script's
zero placeholders at positions 45 and 54 of the partial limit denote missing
entries and are never checked as prescribed data.

Source: https://arxiv.org/abs/2609.10403v1
