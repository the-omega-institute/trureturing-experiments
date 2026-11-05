# Hamzeh–Ashrafi order-supergraph Hamiltonicity

Run with Python 3 and the standard library:

```sh
python3 check.py
```

`output.txt` is the complete deterministic output. All evidence uses integers.

The source is A. Hamzeh and A. R. Ashrafi, “Some remarks on the order supergraph
of the power graph of a finite group”, International Electronic Journal of Algebra
26 (2019), 1–12, DOI [10.24330/ieja.586838](https://doi.org/10.24330/ieja.586838).
Conjecture 2.5 says that the main supergraph of a finite nonabelian simple group
is not Hamiltonian. The simple graph contains the identity; adjacency means
comparability of element orders under divisibility.

The check reproduces Theorem 2.1 for cyclic groups of orders 3 through 32,
Theorem 2.3 on C6, and the obstruction in Theorem 2.2 on A5. Direct enumeration
of A5 independently checks the cycle-type counting routine. The A10 spectrum
also agrees with Theorem 2.31 of the distinct earlier paper
[arXiv:1707.07251v1](https://arxiv.org/abs/1707.07251v1).

For A10, every even integer partition is counted by
`10! / product(j**a_j * a_j!)`, including fixed points. The resulting 13 fibers
sum to 1,814,400. The displayed cyclic order word is expanded into disjoint
intervals of fiber indices. The program checks exact coverage and all 1,814,400
edges, including closure. Any enumeration of each actual order fiber transports
this certificate to the graph on group elements. Permutations themselves are
not exhaustively listed. Simplicity of A10 uses the classical alternating-group
simplicity theorem; the program does not prove simplicity or provide a Lean
proof of the conjecture's negation.

A subset dynamic program checks that the order-label quotient is itself
non-Hamiltonian. The order-21 permutation (1 2 3 4 5 6 7)(8 9 10) is explicitly
checked to be even and adjacent to every order-7 element. Further computations
check the proposed connectedness and spanning-tree capacity criteria for
5 ≤ n ≤ 34, including the five cases n = 6k+4 for 1 ≤ k ≤ 5. These computations
make no unbounded family claim.
