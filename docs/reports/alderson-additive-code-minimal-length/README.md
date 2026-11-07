# Additive code minimal length: exact check

Run `python3 check.py`. Python's standard library suffices; all arithmetic is exact.
The stored output is in `output.txt` and ends with `ALL_CHECKS_OK`.

The checker reproduces the thirty-plane construction of Tim L. Alderson,
*On the Maximality of Additive Codes*, arXiv:2607.22297v1, Section 7,
including its 64 codewords, the distances 24 and 30, and its nonlinear extension.
It independently enumerates all 1395 three-dimensional subspaces of F_2^6.

Replacing the source's sixfold cover by the ten planes from a one-factorization
of K_6 gives a faithful F_2-linear (10,2,8)_(8/2) code. Its nonlinear extension
has parameters (11,2,9)_8. Every three-dimensional subspace meets the set of
minimum-weight differences, which blocks any appended F_2-linear coordinate,
including maps of rank less than three. Puncturing an arbitrary additive
extension would give such a linear coordinate, so the check supports additive
maximality with no faithfulness restriction on the extending coordinate.

An alternating-group orbit gives a threefold cover by fifteen planes. Combining
the twofold and threefold covers supplies the proposed family
(5*lambda,2,4*lambda)_(8/2) for every integer lambda >= 2. The checker verifies
lambda = 2 through 32; the unbounded family requires the algebraic combination
argument in the preregistration. It also checks that no three singular planes
are pairwise disjoint, obstructing a positive onefold cover within this fixed
forbidden-set construction. It does not establish that 10 is globally minimal.

These results are exact computational evidence for a preregistration outline,
not a Lean proof or a formal settlement. The classical regular-cover ingredient
is not claimed new: see Cossidente, Marino, Pavese and Smaldore,
*On regular systems of finite classical polar spaces*, arXiv:2103.09336v1,
Theorem 3.3(2), DOI 10.1016/j.ejc.2021.103439.

The program is original code by this research seat, released under the
repository's license; no third-party implementation is copied.
