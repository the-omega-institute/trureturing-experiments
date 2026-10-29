# Bala congruences for OEIS A193237 and A193236

`check.py` uses Python 3 integer arithmetic and `fractions.Fraction`. It encodes
Hanna's two generating functions literally and computes the constant-one
reciprocal cube roots by the formal differential equation `3 F A′ + F′ A = 0`.
Every coefficient division is exact. A separate binomial expansion and direct
multiplication check `A³ F = 1` through degree 100.

The 21-term prefixes quote OEIS A193237 revision 18 and A193236 revision 20.
Both entries contain Peter Bala's mod-5 and mod-7 conjectures. Data attribution:
Paul D. Hanna, Peter Bala, G. C. Greubel and The OEIS Foundation;
[OEIS End-User License Agreement](https://oeis.org/LICENSE).

Run the standalone check with `python3 check.py`. To reproduce the stored
`output.txt`, download the b-files from
<https://oeis.org/A193237/b193237.txt> and
<https://oeis.org/A193236/b193236.txt> into a temporary directory and run
`python3 check.py --bfile-dir /path/to/temporary/directory`.

Both b-files agree at all indices 0–1000. For each sequence, the exact congruence
checks cover n=0–700 modulo 5 and n=0–500 modulo 7. The finite support and
Frobenius transfer identities are also checked through degree 3505.
The final line is `ALL_CHECKS_OK`. These finite checks do not prove the
unbounded conjectures; Lean proof is outside this preregistration check.
