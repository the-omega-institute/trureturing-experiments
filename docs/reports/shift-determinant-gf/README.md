# Shifted 3 by 3 determinant checks

`check.py` uses Python 3's standard library, integer arithmetic and
`fractions.Fraction`. Run:

```sh
python3 docs/reports/shift-determinant-gf/check.py
```

`output.txt` is the captured successful output. Its final line is `ALL_CHECKS_OK`.

The program reproduces all 45 terms of [OEIS A097399, revision 21](https://oeis.org/A097399),
including the displayed determinant 86 and the value 412 for entries 1 through 9.
It checks the endpoint and slope bounds on all 362880 arrangements, verifies the
normalization on every arrangement, and checks all 5040 representatives using
independent rational Gaussian elimination. It expands each representative's
determinant polynomial exactly, checks the four attaining arrangements, directly
checks shifts 0 through 45, and checks the generating-function numerator through
degree 1000. Finite coefficient checks are not a proof for every degree.

The intended candidate is mathematically consistent with these checks. The
literature screen, however, found Brother U. Alfred, *Note on Third Order
Determinants*, The Fibonacci Quarterly 3(1), 59-60,
[DOI 10.1080/00150517.1965.12431457](https://doi.org/10.1080/00150517.1965.12431457),
[original journal scan](https://www.fq.math.ca/Scanned/3-1/alfred2.pdf).
The body gives the affine shift law, maximal shift coefficient 81, the branches
412+80t and 405+81t, and their crossing at t=7 for entries 1+t through 9+t.
These are precisely the positive-entry part of the candidate after t=n-5.
The last paragraph changes to the first entry m but fails to replace t by m-1;
its printed formulas must not be used literally. The note cites a finite table
and an independent investigation without printing the latter.

This report makes no novelty claim and contains no Lean proof. The source
sequence data are attributed to OEIS under its [end-user license](https://oeis.org/LICENSE).
