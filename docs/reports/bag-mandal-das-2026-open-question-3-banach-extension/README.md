# Bag–Mandal–Das: Banach extension of Open Question 3

Source: P. Bag, P. N. Mandal and P. K. Das, *Asymptotic Expansive and
Hyperbolic Linear Operators*, arXiv:2610.01196v1, Theorem 4.11 and final
Open Question 3. <https://arxiv.org/abs/2610.01196v1>

Run with Python 3.12 or newer, without optimization flags:

```sh
python3 docs/reports/bag-mandal-das-2026-open-question-3-banach-extension/check.py
```

The program uses only the standard library and exact integers and
`fractions.Fraction`. Its output is stored in `check.out`.

- Fidelity: the source example `T = diag(0,2,3)` gives the stated iterate
  and squared Euclidean norm on 1,500 integer-vector and positive-time cases.
- Difference of powers: both orders of the commuting factorization
  `(lambda I - T) sum_k lambda^(N-1-k) T^k = lambda^N I - T^N`
  are checked on 3,906 cases, including the zero-dimensional case.
- Resolvents: all 81 two-by-two matrices with entries in `{-1,0,1}/8`
  have maximum-row-sum norm at most `1/4`. Their powers and all nonsingular
  inverse matrices pass 2,916 contracting and 1,440 expanding rational
  resolvent checks at the listed test parameters.
- Source consistency: for `diag(1/2,2)`, exact squared coordinate weights
  make the input and its image at time `n` both have Euclidean norm one,
  for `1 <= n <= 32`. This concerns a separate converse assertion and is
  not a refutation of the nominated Question 3.

The expected settlement uses a general continuity argument: a homogeneous
norm gap at one common positive iterate cannot mix contracting and expanding
unit vectors. Surjectivity then permits the inverse in the expanding branch.
The finite checks above support identities and examples; they do not prove
the Banach-space theorem, implement its continuity argument, or establish a
complete literature clearance. No Lean proof is included.
