# Hanna's A381361 modulo-five tail

Paul D. Hanna's [OEIS A381361](https://oeis.org/A381361), revision 25, states:
“Conjecture: for n > 4, a(n) == 0 (mod 5).” Its definition is the formal EGF
A = exp(integral abs(A^(-2))), with unsigned coefficients taken coefficientwise.

Run `python3 docs/reports/hanna-a381361/check.py` from the experiments repository.
The program uses only Python integers and fractions. It reproduces all 501
published terms at indices 0 through 500 using the triangular integer recurrences
for A, A squared and its reciprocal. An independent ordinary rational-coefficient
iteration checks the literal source equation through degree 24. The program also
checks the modulo-five tail of all three sequences and the ten small carry pairs
used by the proposed induction. The stored result is `check.out`.

These finite checks establish source fidelity and support the proposed induction;
they do not prove the unbounded conjecture. No Lean proof is included.

`source-b381361.txt` is the unmodified [OEIS b-file](https://oeis.org/A381361/b381361.txt)
by Paul D. Hanna, with SHA-256
`93779794cd57e856400e4156f8eb50cc80971ebed6549b7c2f6647ce6f71dfd4`.
The OEIS data are licensed under
[CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/), as specified by the
[OEIS End-User License Agreement](https://oeis.org/LICENSE).
The independently written program is Apache-2.0, matching this repository.
