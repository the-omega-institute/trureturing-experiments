# Class-two unipotent-action vanishing

For Diego Tuzzolo, arXiv:2610.06909v1, Problem 6.3 (final question).
Run `python3 check.py` with Python 3; only the standard library is required.

The program uses integers modulo the stated prime. It constructs W as the
quotient of the exterior square by K, constructs the Jacobi image in V tensor W,
and computes the action of every basis map V to W on every basis vector of K.
The assertion checked is that the action has rank zero exactly when K is zero
or the Jacobi map is surjective. It also checks that the contractions of K span
a space S whose exterior square contains K.

The exhaustive set comprises all subspaces of the exterior square for
1 <= dim(V) <= 3, over each of the fields with 3, 5 and 7 elements: respectively
31, 67 and 119 cases. A further 240 reproducible cases per prime have dimensions
4 and 5, using seed 261006909. Abelian and zero-dimensional quotient cases are
included as checks of the linear statement's boundary; the proposed group
statement concerns groups of class exactly two.

For source fidelity, Example 5.6 is reproduced at p = 5, 7: dim(K) = 2,
dim(Q) = 2 and rank(Omega) = 1. Proposition 5.7(1) with M = 2 is reproduced
at the same primes: dim(Q) = 4 and rank(Omega) = 4.

Separate exact noncommutative polynomial calculations check the middle tensor
coefficient 2, the inversion signs of cubic displacements from symmetric and
alternating quadratic corrections, and characteristic-three cube polarization.
The output is in `output.txt`; its final line is ALL_CHECKS_OK.

These finite checks neither prove the unrestricted assertion nor enumerate
algebra automorphisms. In particular, they do not establish the
characteristic-three identification of the group-algebra action with Omega.
The proposed positive settlement requires a general argument and its later
Lean formalization. No Lean proof is included here.
