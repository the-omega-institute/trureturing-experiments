# Mathar's two recurrences for OEIS A192744

`check.py` reproduces the five polynomial reductions and all 22 listed terms in
[OEIS A192744](https://oeis.org/A192744), revision 16 (16 December 2015).
The sequence was entered by Clark Kimberling; the conjectures were posted by
R. J. Mathar in 2014 and 2015. Source data and examples are used under the
[OEIS End-User License Agreement](https://oeis.org/LICENSE).

Run with Python 3, without third-party packages:

```sh
python3 docs/reports/mathar-a192744-recurrences/check.py
```

`output.txt` is the complete output. The script uses integer coefficient lists
to build the source polynomials, generic monic long division by `X^2-X-1`, and
the literal iterative substitution defined in
[A192232](https://oeis.org/A192232). These agree through index 200 and agree
with a separate two-coordinate construction. It checks the 2014 recurrence
at all 197 indices from 4 through 200, and the 2015 recurrence at all 198
indices from 3 through 200, together with their factorial-forcing identities.
The conjectured recurrences are never used to generate the test sequence.

The final line is `ALL_CHECKS_OK`; the program exits 0. These are finite exact
checks supporting source fidelity and a preregistration outline. They are
neither an unbounded proof nor a Lean verification of either conjecture.
