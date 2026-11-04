# Miska–Ulas d-cycle Hankel determinant checks

`check.py` uses Python 3's standard library and exact integer/rational arithmetic.
Run from the repository root:

```sh
python3 docs/reports/miska-ulas-hankel-dcycle-vanishing/check.py
```

The source is Piotr Miska and Maciej Ulas, *On some properties of the number of
permutations being products of pairwise disjoint d-cycles*,
[arXiv:1904.03395v1](https://arxiv.org/abs/1904.03395v1), section 8 and
Conjecture 9.5(3); published in *Monatshefte für Mathematik* 192 (2020), 125–183,
[DOI:10.1007/s00605-020-01397-5](https://doi.org/10.1007/s00605-020-01397-5).

The polynomial recurrence is checked against the source's factorial expansion
in 539 cases and its weighted permutation interpretation by direct enumeration
in 63 cases. The conjectured determinant vanishing pattern is checked in 960
cases with 3 ≤ d ≤ 12, 1 ≤ n ≤ 24 and x in {0,1,2,3}. The predicted factorization
has 77 nonzero cases with positive factors and 163 zero cases. Polynomial
constancy is certified for 48 fixed parameter pairs by sufficiently many exact
integer evaluations, using the degree bound n(n−1). Auxiliary checks cover 420
shifted factorial/Pochhammer identities and 405 entries of the positive LDL
identity. Bareiss determinants are cross-checked against Leibniz sums on ten
matrices, including singular matrices and the empty matrix.

`output.txt` is the program's output. Every asserted equality uses integers or
`fractions.Fraction`; no floating point supplies evidence. A successful run
prints `ALL_CHECKS_OK` and exits 0. These finite checks do not establish the
universal conjecture, literature novelty, a Lean proof, or a frozen result.
