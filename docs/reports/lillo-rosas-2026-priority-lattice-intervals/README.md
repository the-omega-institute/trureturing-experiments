# Lillo–Rosas arXiv:2603.28905 §7: principal filters and ideals of the priority lattice

Numerical check for the settlement preregistered in the main repository's preregistration issue for this paper. The kernel-checked results and their dossier are in the main repository; this directory holds only the experimental program.

Run from this directory: `python3 check.py 6` (requires `networkx`).

It enumerates the priority lattice Π(n) for n ≤ 6 (priority forests on {0,…,n} ordered by edge inclusion, plus a top). It counts the principal filters [P, 1̂] isomorphic to some Π(m) with m ≤ n, and the principal ideals [0̂, P] isomorphic to some Π(m) with 1 ≤ m ≤ n and with 0 ≤ m ≤ n.

Files:
- `check.py`: sha256 `b71ca8894f15984e5eca4712f947556f024d3d2ddfc7ebe68142090fc7f562ec`
