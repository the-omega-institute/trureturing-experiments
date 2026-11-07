# Arvind–Hunter-Jones dual-unitary LPM gap

Exact diagnostics for Conjecture 1 of Rahul Arvind and Nicholas Hunter-Jones,
*Finite-bath projected ensembles in dual-unitary circuits*, arXiv:2610.06754v1.
Source: https://arxiv.org/abs/2610.06754v1, Definition 10, Lemma 3, Conjecture 1.

Run with Python 3 and SymPy:

```sh
python3 check.py > output.txt
```

The source-fidelity checks reproduce dual unitarity, traceless one-site operator
spreading, the m=1 LPM gap of 1, matched-permutation invariance, and the dependence
of the six qubit permutation operators at m=3 (their span has dimension 5).

The settlement checks verify the parameter-dependent Pauli conjugations and all
brackets needed for the proposed Lie-generation argument as Laurent rational
identities with t=exp(iJ). They check the Pauli Gram matrix exactly. At m=2,
strict contraction on the 30-dimensional bad sector is checked for t=(3+4i)/5,
t=(5+12i)/13 and a local dressing of the first gate. A Hermitian Schur-complement
factorization with strictly positive rational pivots certifies each contraction.
All arithmetic is exact; no floating point comparison is evidence.

These tests do not establish the unbounded-replica assertion. Its proposed route
uses Lie generation, compact-group bounded exact generation, Schur–Weyl duality
and Haar-projection leakage. No Lean proof of the conjecture or its negation is
included here. The normal form is credited to Bertini–Kos–Prosen,
arXiv:1904.02140, Eqs. (23)–(24); uniform finite generation is credited to
D’Alessandro, arXiv:quant-ph/0111133, Theorem 2. The code is independently written
and does not copy source code from those works.
