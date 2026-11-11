# Kurkov's shifted parity transform

Exact finite checks for the second new conjecture in [MathOverflow question 506710](https://mathoverflow.net/questions/506710/pair-of-identities-for-a125790), using the source's definitions of the signed odd-power moments T and the array P from [OEIS A125790](https://oeis.org/A125790).

Run with Python 3.12 or newer:

```sh
python3 docs/reports/kurkov-mo506710-second-parity-transform/check.py
```

The standard-library program uses integers and rational fractions. It reproduces the first five displayed rows of A125790 and A391910, checks P against an independent direct recurrence table on 189 pairs, and checks the accelerated T moments against their literal binary-sign sums on 341 pairs. Newton polynomials obtained by discrete antidifferentiation of the defining P recurrence avoid exponential table growth. The T expansion splits off one binary digit; neither acceleration assumes a transform identity.

Both previously published transforms pass for 0 ≤ n ≤ 30 and 0 ≤ m ≤ 10 (341 pairs each). The target shifted parity transform passes for 1 ≤ n ≤ 30 and 0 ≤ m ≤ 10 (330 pairs). The stored output is `output.txt`; its final line is `ALL_CHECKS_OK`.

These finite computations establish transcription fidelity and check a bounded range. They do not prove the unbounded conjecture.
