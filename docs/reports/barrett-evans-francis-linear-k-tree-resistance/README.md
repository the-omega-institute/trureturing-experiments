# Straight linear k-tree resistance: exact finite checks

The graph has vertices 1 through n, an edge exactly when 0 < |i-j| ≤ k,
and unit conductance on every edge. The target is the fixed-k endpoint
resistance increment limit 6/[k(k+1)(2k+1)].

Sources: Barrett, Evans and Francis, *Resistance distance in straight linear
2-trees*, arXiv:1712.05883v1, the definition and conjecture numbered 43 and 44
in the arXiv PDF; Evans and Francis, arXiv:2108.07942v2, Conjecture 4.5.
The arXiv v1 PDF numbering is used here.

Run `python3 check.py` with Python 3. No external packages are required.
All evidence uses integers and `fractions.Fraction`. `output.txt` is its
complete standard output. Failed assertions or singular linear systems
produce a nonzero exit; success ends with `ALL_CHECKS_OK`.

The fidelity check reproduces the source's formula
R(n,2) = (n-1)/5 + 4 F(n-1)/(5 L(n-1)) for n=3 through 30.

For each k=1 through 8 and n=2 through 2k+8 (128 cases), the program checks:

- agreement of a grounded-Laplacian solution with the independently assembled
  interval-Gram formula R(n,k) = 1^T Q^(-1) 1;
- the bounds I ≤ Q ≤ D I by exact positive-semidefinite Schur elimination,
  where D = sum(d^2, d=1..k), and the bandwidth bound;
- the exact correction R(n,k) = (n-1)/D + (1^T c + c^T Q^(-1)c)/D^2,
  where c = D 1 - Q 1;
- the explicit, disjoint endpoint defects when n-1 ≥ 2k.

The localization check compares 49 same-endpoint moments at successive
matrix sizes and verifies 98 opposite-endpoint zero moments, for k=2
through 8 and powers 0 through 6. It uses the integer matrix B=D I-Q,
which differs from A=I-Q/D by the scalar factor D.

These computations check finite identities and source fidelity. They do
not prove convergence or any infinite-family theorem. The printed finite
increment errors can have either sign; no monotonicity is asserted.
