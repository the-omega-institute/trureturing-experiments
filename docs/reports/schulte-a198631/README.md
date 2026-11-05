# Schulte's A198631 weighted type-B Stirling identity

Run `python3 docs/reports/schulte-a198631/check.py` from the experiments repository root. Python 3's standard library is sufficient. The output is in `output.txt`.

The program builds A039755 by its defining integer recurrence and logistic exponential-series coefficients by exact rational series inversion. It reproduces the source's A039755 rows through 4 and A198631 numerators and displayed rational coefficients through 12. It checks the conjectured weighted sum, numerator normalization and exact denominator for every index from 0 through 80. An independent truncated power-series product checks the known column generating function through degree 16. All arithmetic used as evidence is integer or `fractions.Fraction` arithmetic.

These checks corroborate finite cases only. They provide neither a proof of the conjecture for all indices nor a Lean certificate.

Sources: Werner Schulte's FORMULA conjecture dated February 16, 2024 in [OEIS A198631](https://oeis.org/A198631), revision 90 (March 12, 2026); definitions in [A039755](https://oeis.org/A039755), [A001147](https://oeis.org/A001147) and [A006519](https://oeis.org/A006519). OEIS data are subject to the [OEIS End-User License Agreement](https://oeis.org/LICENSE).
