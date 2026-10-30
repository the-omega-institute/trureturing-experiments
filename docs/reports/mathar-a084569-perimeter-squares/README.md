# Mathar's A084569 perimeter-square conjecture: exact checks

Run `python3 check.py` with Python 3 and its standard library. The stored output is
`output.txt`; a successful run ends in `ALL_CHECKS_OK`. All arithmetic is integral
or uses `fractions.Fraction`.

The check reproduces the three magic-sum-4 squares in the OEIS comment and the
nine terms printed in Example 7 of Mathar's software README. It independently
enumerates all eight positive labels for source indices 0 through 4, compares
these with corner reconstruction, and enumerates corners, actual C4 orbits and
all rotation fixed sets for indices 0 through 16. It checks every stabilizer
size, reconstruction and equivariance on this range. The source's triple sum
and quartic formula agree for indices 0 through 40. At index 3 the C4 count is
21 and the D4 count is 19, distinguishing rotations from reflections. A one-slot
perimeter shift is also checked to be an invalid substitute for a quarter-turn.

These are bounded checks supporting a proposed general proof. They do not prove
the universally quantified conjecture and are not a Lean probe.

Sources:

- [OEIS A084569](https://oeis.org/A084569), revision 19; COMMENTS conjecture by
  R. J. Mathar, 8 March 2025; FORMULA triple sum and quartic expression.
- R. J. Mathar, [Generating Perimeter-magic Polygons](https://doi.org/10.5281/zenodo.15001366),
  software version 0.069, README of 10 March 2025, Example 7.

The check is independently written; source software and PDFs are not copied into
this report directory.
