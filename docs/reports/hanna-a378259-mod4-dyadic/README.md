# Hanna A378259: exact finite checks

Paul D. Hanna's [OEIS A378259](https://oeis.org/A378259), revision 15
(July 3, 2025; entry created December 4, 2024), conjectures that for n > 1,
a(n) is 2 modulo 4 precisely when n = 2^k + 1, and is 0 modulo 4 otherwise.
The defining generating function is A(x) = sum_{m>=1} A(x^(2m))/A(x^m),
A(0) = 0, A'(0) = 1, with a(n) = 4^(n-1) [x^n] A(x).

Run with Python 3 (standard library only):

```sh
python3 docs/reports/hanna-a378259-mod4-dyadic/check.py
```

`check.py` reconstructs 500 terms from the original rational equation using
`fractions.Fraction`, then independently checks the quotient product and the
coefficientwise Lambert equation. All 500 terms equal the attributed source
table, the 28 displayed source terms agree, and the first negative term is
at n = 224. It checks the proposed scaled square equation through degree 60,
even integrality and the mod-4 conjecture through n = 500, and the mod-8 and
characteristic-two residual identities through degree 499.
`check.stdout` is the complete output; its final line is `ALL_CHECKS_OK`.
These finite checks do not prove the unbounded conjecture. No Lean proof is
included.

`source-b378259.txt` is Paul D. Hanna's table of n, a(n), n = 1..500, obtained
from <https://oeis.org/A378259/b378259.txt>. Source-table SHA-256:
`7478c1c8db8943c9a4ebd9b098070ebe6eff842504602f2c36e6556192405aa6`.
OEIS material is subject to the [OEIS End-User License Agreement](https://oeis.org/LICENSE).
The source identifier, attribution and exact table are retained to make the
fidelity check reproducible without network access.
