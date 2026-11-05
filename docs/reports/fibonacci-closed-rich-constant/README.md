# Fibonacci closed-rich constant: exact finite checks

The target is Conjecture 5.3 in Anuran Maity and Svetlana Puzynina,
*Bounds on the closed-rich constant of infinite words*, arXiv:2605.19535v1.
The source uses f[-1] = b, f[0] = a and f[k] = f[k-1]f[k-2].
Closed factors are distinct words, include the empty word, and may have
overlapping prefix and suffix occurrences of their witnessing border.

Run with Python 3:

```sh
python3 check.py
```

`check.out` records the exact output. The checker uses integers and
`fractions.Fraction`; it does not use floating point or numerical tolerances.
In the quadratic field with sqrt(5), the conjectured constant simplifies
exactly to (2 + 5 sqrt(5))/121. Comparing a rational r against this constant
uses the two rational inequalities 121r - 2 >= 0 and (121r - 2)^2 >= 125.

The check reproduces all 33 rows of Table 1's exact prefix counts and minima.
It then checks the conjectured lower bound on every factor of lengths 1
through 64: 2,144 distinct nonempty factors in total. The least observed ratio
is 91/768, attained at length 48 with 273 closed factors. Two definitions of
closedness, one excluding internal occurrences and the other counting every
overlapping occurrence of each border, agree on all 2,144 words.

Factor coverage follows from the substitution a -> ab, b -> a. Its legal
adjacent pairs are aa, ab and ba: the substitution creates no bb, and all
three pairs occur in f[3] = abaab. For each length n, choose a substitution
level k with both tiles f[k] and f[k-1] at least n long. Every n-window meets
at most two consecutive tiles, so the three substituted pairs cover every
possible factor. All three substituted pairs occur in f[k+3], a prefix of
the infinite source word. The program compares these two finite sets and
checks that their size is n+1. This coverage argument is ordinary
mathematics here; it is not a newly compiled Lean theorem.

For k = 7 through 12, the checker also verifies that
f[k-3] f[k] f[k-3][:-2] is a palindromic source factor with length
F[k] + 2F[k-3] - 2 and exactly 1 + F[k-3]F[k] distinct closed factors,
as in the proof of source Theorem 4.19. At length 48 it verifies the
reported unique minimizer. The family test reaches length 553.

These are fidelity checks and bounded evidence consistent with the
conjecture. They do not prove the universal lower bound, global optimality
of the special factors, or equality of the infinite infimum. The separate
first-difference Conjecture 5.2 is not the target.
