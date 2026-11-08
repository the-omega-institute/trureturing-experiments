# Modal groups: Question 5.5 checks

Source: Wojciech Aleksander Wołoszyn, *Modal group theory*, arXiv:2605.14197v1,
Question 5.5, Theorem 2.4 and Corollary 2.5. The expected answer is negative
for the frame of all groups and embeddings.

Run `python3 check.py`; its deterministic output is in `output.txt`.
The program uses only the Python standard library and exact integer, tuple,
set and Boolean computations.

The fidelity checks reproduce the introduction's nonabelian-product example
and test Theorem 2.4 by actual finite overgroups. For every ordered pair of
generators and every element in eight finite groups, doubled left translation
embeds the group in a symmetric group. Swapping the two copies precisely over
the generated subgroup centralizes the generators and fails to centralize
exactly the elements outside that subgroup. This is a finite witness check of
the membership formula, rather than an assumption about centralizers inside
the original group.

The refutation-route checks verify 44,816 reduced words in the two projected
HNN free bases, 64 solved HNN relations expressing the enumerated elements in
two generators, and the exact modal pattern of the countability button.
For the latter, the mathematical equivalence `diamond Gen2 iff countable`
gives two cardinality classes. Uncountability is false in the countable class
but becomes necessary in the uncountable class. Finite diagonal diagnostics
also check the construction used in Cantor's argument.

These computations do not prove the infinite free-basis statement, the
countable two-generator embedding theorem, Cantor's theorem, or existence of
countable existentially closed groups. Those infinitary bridges, and the
translation from this cardinality pattern to the actual group semantics,
require the mathematical proofs in the preregistration and future Lean work.
No finite group here is claimed to be an existentially closed counterexample.
The checks are supporting evidence, not a kernel-checked settlement.
