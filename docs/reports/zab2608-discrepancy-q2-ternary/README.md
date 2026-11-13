# Ternary discrepancy question 2

`check.py` is a deterministic Python 3 exact-arithmetic check for the proposed
refutation of Question 2 of Zabokritskiy (Yohananov), arXiv:2608.01134v3.
It verifies the Hamming-ball volumes and intersection table, the seven-word
one-error-correcting code, the normalized distance-four-pair clique bound, the
counterexample distance histogram, and the strict discrepancy gap.  It uses
integers and `fractions.Fraction` only; no floating-point arithmetic or
optimizer is used.

Run:

```text
python3 check.py
```

The canonical output is in `output.txt`; a successful run ends with
`ALL_CHECKS_OK` and exits 0.
