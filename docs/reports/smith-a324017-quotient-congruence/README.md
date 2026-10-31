# Smith A324017 quotient congruence

Exact numerical checks of the second conjecture in Davis Smith's
[OEIS A324017](https://oeis.org/A324017/internal), revision 91
(16 February 2025). The entry is dated 28 March 2019; the second conjecture
was added in revision 40 on 29 March 2019, initially with the parameter
names interchanged. The current claim is: for natural numbers n > 1 and
m >= 2, floor(((2*n - 1)^^m)/(2*n)) is congruent to 2*(n - 1) modulo 2*n.
The source defines x^^0 = 1 and x^^(h+1) = x^(x^^h).

Run with Python 3, without third-party packages:

```sh
python3 docs/reports/smith-a324017-quotient-congruence/check.py
```

`check.out` stores the output. All evidence uses integers. The generic
modular tower evaluator factors the modulus into prime powers, reduces
exponents by Euler's theorem only for units, handles nonunits with an
exact valuation cutoff and capped exponent, and combines residues by CRT.
It does not use the conjectured formula. It is cross-checked against
3,584 fully constructed tower/modulus cases and against capped exact towers.

The fidelity checks reproduce 24 displayed array entries and both source
quotient examples. The settlement screen checks 202 literal integer
quotients and all 2,189 pairs with 2 <= n <= 200 and 2 <= m <= 12. It also
checks 39,999 cases of the proposed odd-exponent residue lemma and confirms
that height one fails for 2 <= n <= 20. These finite results do not prove
the assertion for unbounded parameters. No Lean proof is part of this report.

The quoted formula and displayed source data are attributed to Davis Smith
and The OEIS Foundation Inc.; the source is distributed under the
[OEIS End-User License Agreement](https://oeis.org/LICENSE).
