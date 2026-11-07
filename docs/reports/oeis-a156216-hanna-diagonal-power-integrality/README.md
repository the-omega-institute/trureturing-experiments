# Hanna A156216: diagonal ghost powers

Paul D. Hanna's [OEIS A156216, second comment](https://oeis.org/A156216/internal)
(revision #3) conjectures that integrality of
`exp(sum(C(n)*x^n/n, n>=1))` implies integrality after replacing `C(n)` by
`C(n)^n`. All series are formal; ordinary coefficients are meant.

Run `python3 check.py`. The program uses Python 3 standard-library integers
and `fractions.Fraction`, with no floating-point arithmetic or dependencies.
The recorded output is `output.txt`; its final line is `ALL_CHECKS_OK`.

The fidelity check reproduces all 14 displayed A156216 coefficients, using
Lucas numbers as ghosts of `1/(1-x-x^2)`. The universal conjecture is tested on
all 3,125 polynomials with constant term 1 and five coefficients in `-2..2`
through degree 24, and all 625 Euler products whose first four signed
exponents lie in `-2..2` through degree 40. Each case checks exact exponential
recurrence, Möbius ledger reconstruction, the independent Euler product, and
local prime-power congruences. Twelve further signed inputs compare both
exponentials with the defining finite sum `sum(H^j/j!)` through degree 24.
A finite modular check includes negative integers and the prime 2.

These results support the expected settlement and test the implementation;
they do not establish the unbounded conjecture. OEIS supplies the quoted
problem and displayed coefficient data; the check program is an independent
implementation.
