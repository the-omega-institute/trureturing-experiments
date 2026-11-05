# Hadamard-closed polyhedra

Exact checks for the proposed refutation of Johnson–Paparella Question 9.1,
arXiv:2409.07682v4, Section 9. The source's Definition 4.1 uses
`C(S) = {x : S diag(x) S^-1 is real and entrywise nonnegative}` and
`P(S) = {x in C(S) : S diag(x) S^-1 e = e}`.

Run with Python 3 and SymPy 1.14.0:

```sh
python3 check.py > output.txt
```

The program uses integers, rational numbers, and exact symbolic algebra. It
reproduces the source's Walsh `H4` row-cone and row-polytope description. It
checks the rational circle/support identities and common symmetric stochastic
eigenbasis symbolically. For every `m = 10,...,64`, an independent exact convex
hull calculation gives `m+4` edges; `m` nominated exposed edges satisfy every
generator inequality, and all 101,860 products of small-polygon generators have
explicit convex certificates on its four axis vertices. The support normals
homogenize to distinct two-dimensional cone faces. Every generating spectrum
has the stated nonnegative matrix realization. The output ends
`ALL_CHECKS_OK`.

These are finite exact checks and symbolic identities. The general
half-space/facet obstruction, geometric extensions from generators to convex
hulls and cones, and the refutation in Lean remain proof obligations. No Lean
proof of the question or its negation is supplied by this experiment.
