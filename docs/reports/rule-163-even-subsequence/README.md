# Rule 163 even subsequence

Exact finite corroboration for Karl-Heinz Hofmann's OEIS A392003 FORMULA
conjecture, “Conjecture: a(2*n) = A266753(n).” The arithmetic sequence is
`B(m) = (7*2^m - 9 + 5*(-1)^m)/6`, for `m >= 1`. A266753 is the binary
integer read from cells `-n,...,n` of the full integer-lattice Rule 163
single-seed evolution, with the leftmost cell most significant.

Run `python3 check.py`. No third-party packages or network access are needed.
`output.txt` is its output. Every arithmetic operation uses integers or
`fractions.Fraction`.

Two independent implementations evolve the same seed: a sparse set of
exceptions from an evolving uniform background, and a padded dense row whose
ends shrink after each step to preserve complete causal cones. They agree on
92,407 cells, including both moving boundaries and exterior cells, through
iteration 300. The check reproduces the 25 published A266753 terms and the
34 published A392003 terms, and checks the row word `100(10)^(n-1)`, the
closed form, the even-subsequence identity, both recurrences, and coefficients
0 through 300 of Barker's generating-function product.

These are finite checks. They provide no proof for unbounded iteration and
no claim of literature novelty.

Sources: [A392003](https://oeis.org/A392003), revision 30, version timestamp
2026-01-02T08:29:06-05:00; [A266753](https://oeis.org/A266753), revision 49,
version timestamp 2025-02-16T08:33:28-05:00. Numerical source prefixes are
included for fidelity checks; OEIS source material is subject to the
[OEIS End-User License](https://oeis.org/LICENSE).
