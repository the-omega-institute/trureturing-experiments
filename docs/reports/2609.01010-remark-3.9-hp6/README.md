# Exact Hoffman–Pereira certificate check

`check.py` checks Baek–Vargas, arXiv:2609.01010v2, Lemma 3.3,
Theorem 3.4, Example 3.11, and the proposed level-one certificates for
all labeled Hoffman–Pereira matrices of orders one through six.

Run with Python 3.12 or newer:

```sh
python3 check.py
```

The check uses integers only. It represents a certificate multiplied by
four, so the half-integral row-square completion is checked without rounding.
For each negative graph it forces all distance-two entries to +1 and enumerates
all remaining 0/+1 entries. Every triangle-free negative graph is enumerated.
For each matrix the whole cubic identity and nonnegativity of the residual
coefficients are checked. Small zero-neighborhood blocks use component sign
vectors whose Gram matrix is dominated entrywise by the block.

`output.txt` is the program output. These are exact experimental checks;
the planned Lean proof uses a general decomposition and square-completion
argument, with local finite cases only as helper lemmas. No Lean proof is
included in this report.

Source: https://arxiv.org/abs/2609.01010v2

The Python program is original work and uses only the Python standard library.
