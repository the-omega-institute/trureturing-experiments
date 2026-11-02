# Barman–Singh symmetric-group code checks

Source: Subrata Barman and Pushpendra Singh, [On Codes in (Generalized)
Symmetric Groups](https://arxiv.org/abs/2606.21928v2), v2, 8 September 2026,
Example 3.8 and Conjecture 3.10.

Run `python3 check.py` with Python 3. The program uses only the standard library,
integers and `fractions.Fraction`. The complete output is `check-output.txt`.

- The coset-count evaluator agrees with literal multiplication of every pair
  for 39 conjugacy-class/Young-subgroup cases in degrees 1 through 4.
- Source Example 3.8 is reproduced by enumerating all permutations of 9 points:
  the classes of cycle types (5,4) and (5,2,1,1) have 18,144 elements each.
  Their union has coefficient 504 against the Young subgroup of shape (7,1,1),
  and coefficient 144 against the subgroup of shape (6,2,1).
- All 8,888 union/partition cases in Conjecture 3.10 for degrees 5 through 7
  fail to give constant coefficients. Every allowed k and every nonnegative r
  are covered, using the necessary exact equation |X| = r [S_n:S_lambda].
  The moment argument is not used to filter these cases.
- The hook-dominance inequality passes 59,119 cases through degree 25.
  The support polynomial identity and its signs pass 2,499 integer cases.

These are finite checks. They do not prove the unbounded conjecture, establish
Young's rule, or verify any Lean proof. The planned positive settlement uses a
general argument, with these computations serving only as fidelity and error
checks. Multiplicity is retained throughout; a product-set computation would
answer a different question.
