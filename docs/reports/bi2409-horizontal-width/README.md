# Horizontal class-transposition width checks

The target is V. G. Bardakov and A. L. Iskra, *Orders of products of horizontal
class transpositions*, arXiv:2409.13341v1, Question 6.2. The expected answer is
negative. This report supplies exact arithmetic checks, not a Lean settlement.

Run `python3 docs/reports/bi2409-horizontal-width/check.py` from the repository
root. Python 3's standard library suffices. The stored output is `output.txt`;
success means exit 0 with final line `ALL_CHECKS_OK`. No floating point is used.

A word is evaluated as Mathlib permutation multiplication: its rightmost factor
acts first. The source's explicitly printed cycle uses right actions, so that
one fidelity check reverses the written factors. This convention does not
change element orders, word lengths, or the proposed width statement.

The checks reproduce valid displayed source orders 3, 4, 6 and 12, and the
Section 3 cycle (2,4,5,3) with residues 0 and 1 fixed. They also preserve an
unrelated source discrepancy: Corollary 5.5 displays order 2 for the product
of the swaps 0(2) with 1(2), and 0(4) with 2(4). Its exact order is 4.
No step of the width argument uses that displayed equality.

There are 120 generators with moduli 2 through 9. All 14,400 ordered pairs
have product orders in {1,2,3,4,6,12}, realizing every value, and obey the
proposed displacement inclusion and cardinal bound. Every word of length at
most 4 on the generators of moduli 2 through 4 is checked (11,111 words).
A common multiple of the moduli gives a complete finite representation of the
integer permutation; sampling positive points alone is not used to infer an
infinite permutation's order or displacement set.

For block reversals of sizes 2 through 128, the program checks the explicit
horizontal-transposition factorization, involution property, exact displacement
set, and fixed points on every residue in seven consecutive blocks, including
negative blocks. For B from 0 through 7 it checks the whole residue block of
size n = 3^B + 1, including the factorization, and verifies exactly n nonzero
displacements. Every such witness is a nonidentity fixed-point-free involution.

The general outline is: a horizontal generator has at most three displacements;
displacement cardinalities multiply at most under composition; a block reversal
of size n has exactly n displacements. Thus the expected lower bound on its
word length is the least k with n <= 3^k. Finite tests do not prove this
unbounded statement. The explicit factorization supplies only the upper bound
floor(n/2); a matching logarithmic upper bound is not established here.
