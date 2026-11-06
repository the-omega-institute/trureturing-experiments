# Lower and upper mmd-expansiveness

Target: Gabriel Lacerda, *Typical expansive measures for infinite-dimensional maps*, arXiv:2610.09482v1, §2.2, first unnumbered Question.

Run `python3 check.py`. The stored result is `output.txt`; the last line must be `ALL_CHECKS_OK`.
The program uses Python integers and `fractions.Fraction`, with no floating-point evidence or external dependencies.

The fidelity check reproduces the two inclusions of source Lemma 2.1 on exact periodic sequences of periods 1 through 4, using the exact identity that the two-sided metric weights sum to 3.

The proposed family has base q at least 3 and digit support equal to the union of the integer blocks (N_(2r), N_(2r+1)], where N_j = 2^(2^j). The checks cover q = 3 through 9, prefix and tail counts through digit 24, four endpoint pairs through N_7, separated orbit words through length 3, finite-window cover budgets through length 8, and periodic tail inclusions for five positive radii. The output records the precise counts and parameter ranges.

These are finite tests of the proposed general argument. They establish no infinite limit, compactness theorem, homeomorphism theorem, or Lean settlement. The intended argument uses endpoint domination, prefix cylinders, independent coordinate choices, and a metric tail fixed before the time limit. The local step retains a full tail alphabet shift inside every positive-radius Bowen ball at the zero fixed point.
