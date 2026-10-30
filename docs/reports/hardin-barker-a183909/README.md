# Hardin–Barker OEIS A183909

Run `python3 check.py` from this directory. The program uses Python's standard
library, integer arithmetic and `fractions.Fraction`. Its output is in `check.out`.

The source is [OEIS A183909, revision 8](https://oeis.org/A183909/internal),
with R. H. Hardin's [76-term table](https://oeis.org/A183909/b183909.txt).
The arrangement definition and data are attributed to R. H. Hardin; the
generating function and recurrence are Colin Barker's conjectures. OEIS source
material is licensed under [CC BY-SA 4.0](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).
The occurrence interpretation and capped-state method are credited to Robert
Israel's [parent-family note](https://oeis.org/A183912/a183912.pdf).

For each represented occurrence, the two summands must occupy distinct other
positions; their values may coincide. Two independent predicates are compared:
multiplicity witnesses, and size-two subset-sum dynamic programming on actual
occurrences after deleting the represented occurrence.

The exact checks cover:

- all 16,384 multiplicity vectors capped at 3, including agreement of both
  predicates, 13,234 valid vectors and the complete finite numerator identity;
- all 50,352 nondecreasing arrangements at source indices 1 through 10,
  reproducing the unique all-zero arrangement at index 2 and the displayed counts;
- the 34 displayed terms and all 76 b-file terms, compared with exact generating
  function coefficients and independent capped-fiber coefficients;
- Hardin's empirical polynomial through index 76 and Barker's recurrence through
  index 76, with nonzero recurrence residual 4 at index 10.

These are experimental results. An all-length conclusion additionally requires
the mathematical cap-invariance, occurrence-to-multiplicity equivalence and
fiber-generating-series arguments. No Lean proof is included.
