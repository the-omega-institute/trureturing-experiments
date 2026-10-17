# OEIS A380680: the mod 4 conjecture is misstated

A380680 is the series Σ_{n∈ℤ} x^(n²)(2 − x^n)^n/(1 − 2x^n)^n. Its conjecture says a(n) ≡ 2 (mod 4) iff n is a square, and 4 | a(n) otherwise. The same wording appears in A380060. Here it fails at once: a(1) = 4, and a(2) = 6.

**Corrected statement.** a(m) ≡ 2 (mod 4) iff m is twice a square, and a(m) ≡ 0 otherwise. The reason: a(m) = 2c(m), and modulo 2 each term x^(n²)(2 − x^n)^n(1 − 2x^n)^(−n) reduces to x^(2n²).

`check.py` (Python 3 with sympy) expands the series exactly. It then checks:
- the published data;
- the two counterexamples;
- the corrected statement for m ≤ M.

Run `python3 check.py 120`. Exit 0 and `RESULT PASS` mean every check holds.
