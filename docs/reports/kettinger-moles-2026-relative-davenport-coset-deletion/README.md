# Relative Davenport constants after coset deletion

Source: Jared Kettinger and Grant Moles, *Elasticity of Orders from the
S-relative Davenport Constant: an Arithmetic Application of a Number-Theoretic
Investigation*, [arXiv:2605.17595v2](https://arxiv.org/abs/2605.17595v2),
Definitions 2.1 and 2.4, Theorem 2.13, Corollary 2.15, Conjecture 2.17.

Run `python3 check.py`. The stored output is `check.output.txt`; success exits
zero and ends with `ALL_CHECKS_OK`. Python's standard library suffices. Every
calculation uses integers and finite sets; no floating point is used.

An S-sum sequence has its total in S; its terms may lie anywhere in the group.
Zero-sum-free excludes all nonempty positional subsequences, including the full
sequence. Repetitions are allowed. The empty sequence contributes only length
zero and does not change the source's zero fallback.

The checker enumerates all zero-sum-free multisets in cyclic groups C_n for
1 <= n <= 12. Sorting loses no information in an abelian group. A branch is
pruned only after it contains a zero-sum subsequence, and every possible next
term is tested even at length n-1. The checker reproduces the source's
singleton evaluations and Theorem 2.13 for every cyclic coset in that range.
For n <= 8 it also checks Lemma 2.5 and Theorem 2.6 for every nonempty S.

A separate literal powerset implementation tests all 341 ordered C_4 sequences
of lengths zero through four. With H1={0,2}, H2=C_4 and g=1, it obtains
alpha={1,3}, beta=C_4, d_beta=3 and d_(beta minus alpha)=2. The same gap is
checked for n=4,6,8,10,12: d_C_n=n-1 and d_(2C_n)=n-2. These finite checks do
not prove the infinite family.

The Corollary 2.15 proof's asserted gcd preservation fails at p=2, a=2, b=0:
adding p^b to g1=1 gives g1prime=2, with gcd(2,1)=1 but gcd(2,4)=2.
This is a substantive failure to preserve extremal totals. The unshifted
subgroup form alpha=H1, beta=H2 survives for cyclic ambient groups; the
checker verifies every such nested pair in its range. No conclusion about
all downstream number-field elasticity formulas follows from this check.

Skałba's classical singleton evaluation is prior art: M. Skałba,
*On numbers with a unique representation by a binary quadratic form*,
Acta Arithmetica 64 (1993), 59–68, Theorem 2,
[doi:10.4064/aa-64-1-59-68](https://doi.org/10.4064/aa-64-1-59-68).
The proposed result is the explicit nested-coset refutation, not a new
calculation of classical cyclic relative Davenport constants.
