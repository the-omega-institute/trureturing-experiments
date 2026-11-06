# Exact checks for twin-free theta separation

`check.py` checks Alfuraidan, Khan and Sauras-Altuzarra, *Metric dimension and maximum eigenvalue multiplicity of graph Laplacians*, arXiv:2610.26875v1, Problem 30.

Run with Python 3 and SymPy 1.14.0:

```sh
python3 docs/reports/aks-problem30-twin-free-theta/check.py
```

All graph matrices are constructed from edge lists. Breadth-first search uses integer distances. Characteristic polynomials, squarefree decompositions, polynomial gcds, symmetry restrictions and kernel computations use exact integer or rational arithmetic. Maximum multiplicity is taken over all eigenvalues; real symmetry identifies algebraic and geometric multiplicity. Both adjacent and nonadjacent twins are excluded.

The fidelity checks reproduce rows 1 and 4 of the source's table labelled `TwinFreeExamples`, with `(dim, mm_L, mm_Q)` equal to `(2,2,1)` and `(2,1,2)`. They also reproduce Theorem 25 for `Theta_{2,2,r}` at `r = 3,5,7,9`.

For `Theta_{3,3,6k+4}`, every `k = 1,...,8` passes the twin test, has metric dimension 2, and has maximum multiplicities `mm_L = 2`, `mm_Q = 1`. Its characteristic-polynomial gcds are `(x-1)(x-3)` for L and 1 for Q. The stated distance formulas and both invariant restrictions pass independently against the full edge-list matrices. The symmetric restriction has simple spectrum; its multiplicities at 1 and 3 are `(1,1)` for L and `(0,0)` for Q.

The program also checks the two symbolic obstruction factorizations, their quadratic remainders and irreducibility, the singular transfer periods, and continuant identities through index 80. `output.txt` is the complete output and ends in `ALL_CHECKS_OK`.

These are exact finite checks and symbolic algebra checks. They support a proposed uniform argument; they do not certify the infinite-family statement, constitute a Lean proof, or establish literature novelty. The program contains no downloaded third-party code.
