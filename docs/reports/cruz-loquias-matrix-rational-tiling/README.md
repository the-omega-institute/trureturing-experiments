# Cruz–Loquias fixed-root tiling check

Run `python3 check.py` with Python 3; only the standard library is required. All arithmetic is integral or rational. The final line is `ALL_CHECKS_OK`.

The fidelity check reproduces arXiv:2607.08732v1 Section 4.2: with its displayed matrices, path 443 reaches (−5,3). The family checks cover odd m from 9 through 49: exact child sets through depth 3, endpoint recurrences through depth 12, dyadic quotient orbits for denominator exponents at most 6 and numerator absolute value at most 50, and integer descent for absolute value at most 300. The four endpoint parity cases are checked for odd m through 101.

The program checks rational width and overlap bounds, including the power-zero scaling printed in equations (16)–(17). Its finite-prefix tail certificates rely on the proposed interval-limit lemma; the program does not establish that lemma or prove the infinite family. The accompanying preregistration supplies the mathematical argument. `check.out` is the program output.

Source: Anjelo Gabriel R. Cruz and Manuel Joseph C. Loquias, *A Matrix Analogue of Rational Number Systems*, arXiv:2607.08732v1, https://arxiv.org/html/2607.08732v1 . V2 removes the expansion-tree section. This report concerns v1.
