# Original Favorov measure: exact finite checks

`check.py` uses Python 3 integers and `fractions.Fraction` throughout. Run it with `python3 check.py`; the recorded output is `check.out`.

The fidelity checks reproduce the time and Fourier support gaps of Favorov's Lemma 3 at alpha = 1/8 for M = 4, 8, 16, 32. They construct nonzero coefficients in the exact cyclotomic quotient Z[z]/(z^(M^2/2)+1). This is a concrete instance of the source's periodic-comb lemma; it does not certify every large M or the complete original measure. The coefficients can be scaled to modulus at most one. The source's additional normalization and isolation choices are not certified by this finite check.

The settlement checks cover 528 rational lattice offsets with rational bounds for both infinite tails, 1,932 exact divided differences, finite rational dipoles with all M <= 8 in the displayed grids, and 64 terms of an explicitly admissible decay sequence. The kernel 1/(1+x^2) is an auxiliary C1 envelope test and is not claimed to be a Schwartz function. The Fourier counting inequality is checked for integer M,R <= 16. All arithmetic and the polynomial exponent identities are exact.

These are consistency checks for the proposed uniform-translate argument. They are not a proof of the infinite-family classification, of a Schwartz convolution statement, or of the existence of all construction choices. No Lean proof is included.

Sources: S. Yu. Favorov, arXiv:2401.01121v1, Lemmas 3–4 and the proof of the unnumbered theorem; J. Mazáč, C. Richard and N. Strungaru, arXiv:2605.23884v1, Section 4.1 after Theorem 4.1.
