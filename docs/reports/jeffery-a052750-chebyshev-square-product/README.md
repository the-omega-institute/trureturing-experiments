# Jeffery A052750 Chebyshev square-product checks

The source is L. Edson Jeffery's conjecture in the COMMENTS of
[OEIS A052750](https://oeis.org/A052750), alongside the explicit sequence
and Wolfdieter Lang's polynomial-discriminant example.

Run with Python 3 and SymPy:

```sh
python3 docs/reports/jeffery-a052750-chebyshev-square-product/check.py
```

`check.out` contains the output. Successful execution ends in `ALL_CHECKS_OK`.

The fidelity checks reproduce the b-file entries at indices 0 through 12,
and the source's explicit degree-four polynomial and discriminant 729.
The conjectured product is checked for every integer n from 1 through 20,
including composite values and repeated prime factors of 2n+1.

The first calculation constructs the source nodes as t^(2k−1)+t^(−2k+1)
in the exact field Q[t]/Phi_(4n+2), identifying t with exp(pi i/(2n+1)).
It evaluates the literal recurrence at every node, checks each proposed
factor identity (4−x_k²) sum_m S_(2m)(x_k)² = 2n+1, and multiplies the factors.
The second calculation uses the resultant of S_n−S_(n−1) and the square-sum
polynomial. All arithmetic is in integer or rational polynomial rings;
there are no numerical roots or floating point comparisons.

These finite checks test transcription and the proposed route. They do not
prove the conjecture for all positive n, and contain no Lean proof.
