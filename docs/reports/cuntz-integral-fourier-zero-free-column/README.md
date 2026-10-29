# Cuntz integral Fourier matrices: zero-free columns

Source: Michael Cuntz, *Integral modular data and congruences*, arXiv:math/0611233v2 (3 June 2008), Journal of Algebraic Combinatorics 29 (2009), 357–387, DOI 10.1007/s10801-008-0139-y. Section 7(5) asks whether a zero-free column can contain an integer other than ±1.

Run with Python 3, using only the standard library:

```sh
python3 docs/reports/cuntz-integral-fourier-zero-free-column/check.py
```

The stored output is `output.txt`. All arithmetic used as evidence is integral or rational (`fractions.Fraction`); square roots in the source example use `math.isqrt` followed by an exact square check. The final line is `ALL_CHECKS_OK`.

The fidelity check transcribes the 13-dimensional congruence multiplication table in Section 4, the example immediately after Definition 4.1. It verifies its defining conditions (1)–(6), all 2,197 integral structure constants, and the author's divisibility counterexample: at source indices (6,6), the symmetric scaled entry is 24 and the degree product is 64, so the former does not divide the latter. The row norms are 1600,1600,1600,1600,400,25,25,25,16,16,4,4,4.

The settlement checks enumerate all square integer matrices of sizes 1–4 with entries in [−2,2], first column one and pairwise orthogonal rows, up to row permutation. Nonzero orthogonal rows imply invertibility. There are respectively 1, 1, 8 and 26 such matrices. For each matrix the program checks every entry of the weighted column Gram identity, the reciprocal-norm sum, the zero-excess balance, and every zero-free column's unit entries. It also checks those identities on the source's 13-dimensional example. This covers a larger class than integral Fourier matrices and does not impose symmetrizability or fusion integrality on the enumeration.

Three exact counterexamples show the limits of dropping integrality, the first-column normalization, or squareness. These finite checks support the proposed negative answer; they do not prove the assertion for arbitrary size or unbounded entries. The source's explicit Fourier definition is Definition 4.1(1)–(4), without a T-matrix. Its diagonal equation determines row norms; the reversed-index explanatory sum printed in Definition 3.1 is not used as an extra axiom.
