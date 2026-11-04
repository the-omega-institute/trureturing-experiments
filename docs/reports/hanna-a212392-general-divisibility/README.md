# Hanna self-composition divisibility: exact finite checks

Run with Python 3 and its standard library:

```sh
python3 docs/reports/hanna-a212392-general-divisibility/check.py
```

The program uses sparse ordinary polynomials with integer coefficients and
`fractions.Fraction` for the primitive. It iterates
f = (x + f(f(x)))^b and independently iterates the reduced-coordinate equation
from Lemma 7.1 of Abbondandolo–Majer, *Gibbs measures on Cayley trees and
homoclinic orbits*, [arXiv:2610.33287v1](https://arxiv.org/abs/2610.33287v1).
Both constructions agree for the first 16 supported coefficients for every
b from 2 through 8.

Fidelity checks reproduce all 12 displayed b=2 coefficients and all 11 displayed
b=3 coefficients in Remark 7.2. The conjecture is the general COMMENTS assertion
in [OEIS A212392](https://oeis.org/A212392) and
[OEIS A212425](https://oeis.org/A212425), attributed to Paul D. Hanna.
The source paper identifies the sequences; it does not state that conjecture.

For each checked b, write q=b^2-1 and a(n)=[x^(q(n-1)+b)]f. The program checks
112 exact divisibilities ((b-1)(n-1)+1) | a(n), with n=1 through 16.
It constructs the rational primitive P(0)=0, P'=(b+1)f, checks that its retained
coefficients are integers, and verifies
P-P(f)=(b+1)xf-b(x+f(f))^(b+1) at every degree through b+15q+1.

`output.txt` is the program's complete output. Its final line is
`ALL_CHECKS_OK`; the process exits 0. These are finite checks supporting a
preregistration outline. They prove neither the unbounded conjecture nor any
Lean statement. No floating point arithmetic or external packages are used.
