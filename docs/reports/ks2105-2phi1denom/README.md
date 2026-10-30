# Kim–Stanton contiguous ratio denominators

Source: Jang Soo Kim and Dennis Stanton, *Three families of q-Lommel polynomials*, arXiv:2105.10096v3, Conjectures 8.3 and 8.4, and its Heine contiguous relation and two continued fractions.

Run `python3 docs/reports/ks2105-2phi1denom/check.py` from the experiments repository root. Requirements: Python 3 and SymPy (the stored output records their versions). Arithmetic uses exact multivariate rational functions over QQ and Python `Fraction`; no floating point is used.

`output.txt` records these finite results:

- The published Heine contiguous relation agrees symbolically in degrees 0 through 4.
- For each conjecture, every denominator-cleared quotient coefficient in degrees 0 through 4 is in Z[a,b,c,q]. Polynomial term counts are respectively `[1, 6, 34, 209, 678]` and `[1, 4, 24, 162, 530]`.
- On five nonsingular exact rational parameter tuples, quotient coefficients agree with independently evaluated q-Nörlund and Heine continued fractions through degree 12. The symbolic polynomials also agree with those specializations through degree 4.

These finite checks do not prove either universal conjecture, minimality of the proposed denominator, positivity, or sharpness of any exponent. The program is original code under this repository’s Apache-2.0 license; it contains formulas and source references, with no bundled third-party source text.
