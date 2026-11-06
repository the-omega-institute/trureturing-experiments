# Hanna variable-power exponential integrality: exact finite checks

`python3 check.py` uses Python 3 integers and `fractions.Fraction` only. It reproduces the NAME generating functions, the printed initial coefficients and logarithmic coefficients of OEIS A181070, A181074 and A181080. The exponential is evaluated by its exact differential recurrence and independently by a truncated Taylor expansion through degree 18.

The supplied b-file tables agree in every row (A181070: 0–120; A181074: 0–75; A181080: 0–90). All three exponential coefficient lists and their integer logarithmic coefficients are integral through degree 140. The program checks every Dold pair p*m <= 140 and 26,240 variable-power congruences for the binomial and negative-binomial coefficient families (n=1–40, k=0–40, primes 2,3,5,7,11,13,17,19). It checks the A181080 reindexing j=n-k exactly. Without the extra +1 in the exponent, the first noninteger coefficient is at degree 6: 45/2, 1481/2 and 327/2 respectively.

`output.txt` is the deterministic output; it ends with `ALL_CHECKS_OK`. These finite checks support a proposed proof route; they do not prove any all-degree conjecture. No Lean proof is included.

Source: Paul D. Hanna, OEIS [A181070](https://oeis.org/A181070) (revision 11), [A181074](https://oeis.org/A181074) (revision 9), [A181080](https://oeis.org/A181080) (revision 6). The adjacent b-files are G. C. Greubel’s tables, retrieved from their corresponding OEIS URLs. OEIS source content is governed by the [OEIS End-User License Agreement](https://oeis.org/LICENSE); fixture hashes appear in the output. The program has no third-party Python dependencies and requires no network access.
