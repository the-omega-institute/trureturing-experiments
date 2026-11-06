# Almost translation-invariant measures on the integers

Source: Zeinab Ashtab and David Fernández-Bretón, *Finitely additive measures on ℤ and additive combinatorics*, [arXiv:2607.26522v1](https://arxiv.org/abs/2607.26522v1), Definition 4.1 and Question 5.2(1).

Run `python3 docs/reports/ashtab-fernandez-breton-2607-26522-q5-2-1/check.py` from the experiments repository root. Python 3's integers and `fractions.Fraction` supply every numerical calculation. The program exits 0 and ends with `ALL_CHECKS_OK`; its complete output is in [check.out](check.out).

The fidelity check reproduces the source's parity calculation immediately following Definition 4.1: an invariant probability assigns the even and odd classes mass 1/2 each. The construction checks condition uniform probability on the subgroup qℤ/Mℤ, for q = 2,…,6, 2 ≤ M/q ≤ 6 and M ≤ 12. All 19,600 subsets and 226,624 shifts in these 12 models satisfy the support, subgroup-invariance and almost-invariance identities, and each model has both required separations. These models describe finite algebras of periodic subsets of ℤ.

The program also checks finite additivity on 3,279 disjoint sample pairs, the symmetric-difference estimate on all 3,590 tested window/subset/shift combinations, and the family values on qℤ, qℤ−1 and 2qℤ for sample lengths 1,…,64. These are finite checks of the general construction. They do not establish existence of the ultrafilter limit or almost invariance on the full powerset of ℤ; those require a general compactness and finite-additivity argument.

Definition 2.3(1)(b) prints an intersection in its disjoint-additivity formula. These checks use the conventional disjoint-union law used in the paper's computations.
