# Bala A370261: exact finite checks

Run `python3 docs/reports/bala370261/check.py` with Python 3.12 or newer. It uses only the standard library, integer arithmetic, rational arithmetic and integer square roots. The saved output is `check.out`; success ends with `ALL_CHECKS_OK`.

The program reproduces all 350 published A370259 b-file values and all 300 A370261 b-file values. It independently evaluates the first-kind recurrence by binary matrix powers, checks the Chebyshev doubling and Pell identities, checks the polynomial congruence modulo 2n and odd integer witnesses for n=1..300, and checks the normalized radicand as an exact rational square. The by-product n divides U_(n-1)(1+c*n) is checked for n=1..100 and every c=-7..7, including odd c. A separate invariant grid covers x=-4..4 and m=0..20.

These are finite fidelity and identity checks supporting a preregistration. They do not prove the universal conjecture and are not Lean kernel verification. No Lean proof is included.

Inputs and attribution:

- `source-fields.json`: mathematical fields from [OEIS A370261 JSON](https://oeis.org/search?q=id:A370261&fmt=json), revision 15, published 24 July 2024, and [OEIS A370259 JSON](https://oeis.org/search?q=id:A370259&fmt=json), revision 17, published 14 March 2024. Both entries were authored by Peter Bala, submitted 11 March 2024. The A370261 offset is `1,2`.
- `b370261.txt`: Paolo Xausa, [Table of n, a(n) for n=1..300](https://oeis.org/A370261/b370261.txt).
- `b370259.txt`: Paolo Xausa, [Table of n, a(n) for n=1..350](https://oeis.org/A370259/b370259.txt).

The input data remain under the OEIS [Creative Commons Attribution Share-Alike 4.0 licence](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement). The check program is original to this report and follows this repository's Apache-2.0 licence. It does not copy OEIS program fields or any third-party generator. Input SHA256 values are printed in the output.
