# Ideal zero-product probabilities and Loewy profiles

Source: Sukrit Chakraborty and Sourav Kanti Patra, *Ideal Zero-Product Probability in Finite Commutative Rings*, arXiv:2608.14618v1, Conjecture 4.11.

Run with Python 3 and its standard library:

```sh
python3 check.py
```

`check.out` contains the complete deterministic output. All arithmetic uses integers or `fractions.Fraction`; the program contains no floating-point calculations.

The checks reproduce the source's values 7/18 and 25/64 for the profiles (2,5) and (3,3) at k=2. Ideal products are also computed from residue sets in Z/(2^n) for 1≤n≤8 and 2≤k≤12 (88 cases). Completeness of the ideal list uses the elementary fact that additive subgroups of a cyclic group are cyclic. A separate exhaustive exponent enumeration verifies the nonzero-product count binomial(k+n−1,n−1) and its slack-coordinate bijection in 30 cases, 1≤n≤6 and 2≤k≤6.

The proposed ordered-scale identity is checked in 420 cases. Every profile of at most five entries from lengths 1 through 8 is included: all 1,287 exact probability signatures at k=2 through 16 are distinct. There are 9,009 exact telescoping checks and 9,002 cancellation checks. Fields (length 1), repeated factors and empty profiles are included.

These are finite checks of the source transcription and proposed argument. They do not establish the infinite identification theorem, the asymptotic limit, or a Lean settlement. The intended positive settlement uses a general minimum-extraction and cancellation argument; finite enumeration supplies no admission basis.

SHA-256 of `check.py`: `bdfefd34e1d580c889be46694bbd1a301228d24ff42a77f131d6c47891d9a93b`.
