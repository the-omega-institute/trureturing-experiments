# Browne's A164660 floor product

`check.py` checks the complete final FORMULA conjecture by Anthony Browne in
[OEIS A164660, revision 19](https://oeis.org/search?q=id:A164660&fmt=text):

\[
\operatorname{num}\!\left(\int_0^1 T_{n+1}(x)\,dx\right)
=\prod_{k=1}^{n}\left(1-2\left(\left\lfloor\frac{k^n}{n}\right\rfloor
-\left\lfloor\frac{k^n-1}{n}\right\rfloor\right)\right)
=(-1)^{n/\operatorname{rad}(n)},\qquad n\ge1.
\]

The source introduces the integral numerators through row sums of the rational
coefficient triangle A164658/A164659. `Fraction` uses the reduced signed numerator
and positive denominator. The radical has `rad(1)=1`.

Run with Python 3; no external packages or network access are needed:

```sh
python3 docs/reports/browne-a164660-floor/check.py
```

`output.txt` is the resulting standard output. The program reproduces the eleven
rational values in A164660's `%e` line and the first five rational triangle rows
in A164658/A164659. It constructs `T_0,...,T_301` by the integer polynomial
recurrence and independently checks every coefficient against the binomial
formula. It then checks both conjectured equalities for every `1 <= n <= 300`,
using literal integer powers and floor divisions, and verifies the floor
indicator, exact count, integral formula, and parity bridge in those ranges.
All 45,150 pairs `1 <= k <= n <= 300` pass. These finite checks do not prove the
unbounded assertion.

The linked [Lang table](https://oeis.org/A164660/a164660.txt) contains a sign typo
in its displayed rational row 7: `56/3` should be `-56/3`. Its numerator table
has `-56`; the defining recurrence gives the coefficient `-112` of `X^5` in
`T_7`, hence the primitive coefficient `-56/3` of `X^6`. The source `%e` row sums
used for fidelity are consistent with the recurrence.

The classical count of nilpotents, `n/rad(n)`, appears in Steven Finch,
*Idempotents and Nilpotents Modulo n*, arXiv:math/0605019v2, Section 2.
That ingredient and the classical Chebyshev antiderivative are prior art;
this program claims no novelty for them and contains no Lean proof.
