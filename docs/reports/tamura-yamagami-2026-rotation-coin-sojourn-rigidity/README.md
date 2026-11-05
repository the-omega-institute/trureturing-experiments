# Rotation-coin positive-sojourn rigidity

`check.py` checks the grouped-path normalization and outgoing-step-at-zero convention of Shunya Tamura and Tomoki Yamagami, *Hadamard Rigidity of Positive Sojourn Time Distributions for Rotation Coins*, arXiv:2609.05033v1, Sections 2--5. It uses Python 3 and SymPy, with exact integers, rationals, complex rationals and polynomial algebra throughout.

Run `python3 check.py`. The stored output is `check.out`; successful execution exits 0 and ends with `ALL_CHECKS_OK`.

The source fidelity check reproduces Section 4's Hadamard time-8 matrices, weights `1/128`, normalizer `3/128` and conditional masses `1/3`. Independent whole-word enumeration and position/sojourn recurrence agree in all nine time-8 bins, using all 70 returning words. The polynomial Gram obstruction is checked identically in real indeterminates `c,s`. Balanced coefficients and their state-independent Gram weights are checked for `m=1,...,8`; exact phase covariance is checked at every even time through 12. The two degenerate coin types and the excursion generating-function inverse identity are also checked.

The small-time checks are finite evidence. The symbolic Gram and generating-function identities are exact algebraic checks. Neither is Lean kernel verification, a proof of the universal classification, nor evidence of literature novelty. Conditional probabilities here normalize separately grouped return-path amplitudes; they do not implement a different measurement-based sojourn law.
