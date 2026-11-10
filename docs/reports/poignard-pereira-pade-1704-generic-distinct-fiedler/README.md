# Poignard--Pereira--Pade generic distinct Fiedler coordinates

`check.py` is an independent exact-arithmetic check for arXiv:1704.01677v1. It builds the weighted three-vertex path with edge weights `(1,2)`, verifies the Laplacian and its characteristic polynomial, identifies the exact Fiedler eigenvalue `3 - sqrt(3)`, and checks that the Fiedler vector and every nonconstant eigenvector have nonzero pairwise-distinct coordinates. It also computes the symbolic discriminant of the nonzero characteristic factor for the two-weight path family.

Run with:

```sh
python3 check.py
```

The final line must be `ALL_CHECKS_OK`.
