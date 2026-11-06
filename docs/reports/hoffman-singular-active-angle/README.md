# Singular active LASSO intersection angle

Independent exact checks for the first future-work item of Faris Chaudhry,
Anthea Monod and Keisuke Yano, *Explicit Bounds on Hoffman Constants for Sparse
Polyhedral Regularizers*, arXiv:2610.15241v1, Section 8. The distance aggregation
is the maximum in Equation (20), rather than a sum of squared distances.

Run with Python 3 and SymPy 1.14.0:

```sh
python3 check.py
```

`check-output.txt` contains the output. The check uses rational numbers and
exact real algebraic expressions throughout. It exits zero and prints
`ALL_CHECKS_OK` after the following checks.

- Source fidelity: Proposition 3, with a nonsingular two-coordinate LASSO KKT
  system. Its affine solution is feasible, and the exact least constant is one;
  49 rational points satisfy the source's maximum-distance identity.
- The singular two-observation family has design rows `(1,1,1)` and
  `(1+e,0,e)`, active signs all positive, regularization parameter one,
  feasible point `(0,1,1)`, and response `A p + (2,0)`. Its null direction is
  `(e,1,-1-e)`. Symbolic identities verify the KKT equation and the determinant
  and trace bounds giving a uniform positive lower bound on the two nonzero
  design singular values for `0 < e <= 1`.
- At `e = 1, 1/2, 1/3, 1/5, 1/10, 1/100`, the least admissible spectral value
  equals `e^2/(e^2+1+(1+e)^2)`. Exact attaining points give the proposed constant
  squared, `2/(1-sqrt(1-eta))`.
- Sixteen additional cells include twelve deterministic random line cells,
  three cells with a two-dimensional active nullspace, and a singleton cell
  with a two-dimensional nullspace. The seed is 261015241. Half the cells add
  an inactive zero design column with the source's subgradient box. All 576
  rational sample points satisfy the proposed bound. Every cell also has an
  exact attaining witness.

For each cell, the spectral calculation enumerates coordinate masks and tests
an eigenvector's signs and its normal-cone inequalities against every vertex
of the feasible polytope. A separate nearest-point calculation projects onto
each feasible affine face and takes the least feasible squared distance.
The selected nontrivial eigenspaces are simple; the program explicitly rejects
a repeated nontrivial eigenvalue rather than silently treating one eigenbasis
as an exhaustive sign-feasibility test.

The proposed universal identity is not established by these finite checks.
The general proof requires reducing the minimum over normal directions to
ordinary eigenvectors of the compressed coordinate projectors. No Lean proof
or negation of the target is supplied here.
