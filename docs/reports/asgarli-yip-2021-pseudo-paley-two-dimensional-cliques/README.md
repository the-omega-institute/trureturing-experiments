# Asgarli–Yip arXiv:2110.07176 Conjecture 5.9: two-dimensional cliques in PP(p⁴, p+1, I)

Numerical check for the settlement preregistered in the main repository. The kernel-checked result and its dossier are in the main repository; this directory holds only the experimental program.

Run from this directory: `python3 check.py 3 5 7` (standard library only).

For each odd prime p, it builds F_{p⁴} with a primitive element g, then enumerates the p² + p + 1 two-dimensional F_p-subspaces V containing 1. A subspace V is a clique in PP(p⁴, p+1, I) for some I with |I| = (p+1)/2 exactly when V∖{0} meets at most (p+1)/2 of the (p+1)-th cyclotomic classes. The program compares that set with the subspaces F_p + aF_p, a = g^{(p+1)k}, k odd.

Files:
- `check.py`: sha256 `4c1c5a9ef3fa872513d822230bde1adac20d87a0427ba4aa5bf1d95b420ab56e`
