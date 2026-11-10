# OEIS A341379, conjecture (2)

Paul D. Hanna's [OEIS A341379](https://oeis.org/A341379), revision 9, defines the
normalized formal EGF by

\[
 A(x)^2=\sum_{r\ge0}\frac{x^r}{r!}\exp(xA(x)^r),\qquad A(0)=1.
\]

Conjecture (2) says that its labelled coefficients are never divisible by three.
These finite exact checks support the proposed stronger residue formula: the
coefficient is two modulo three at indices congruent to five modulo six, and one
otherwise. They are not an all-index proof or a Lean kernel verification.

Run with Python 3.12 or newer:

```sh
python3 check.py
```

The default command downloads Kotesovec's [b-file](https://oeis.org/A341379/b341379.txt).
For a supplied identical copy, use `python3 check.py --bfile /path/to/b341379.txt`.
The b-file content hash is recorded in `check.out`; no b-file content is bundled.
The source OEIS data retain their [CC BY-SA 4.0 terms](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).
The check program is independently authored for this report.

The program uses integer binomial convolutions and rational ordinary-series
arithmetic. It reproduces the displayed coefficients and squared-series table,
checks the literal exponential equation through degree 24, and compares all 126
published b-file residues. An unreduced-power source recurrence through degree
64 cross-checks a separate Frobenius-based recurrence through degree 2000. The
large recurrence does not take the proposed residue formula as input. It also
checks the nonlinear reduction and square-zero correction on that finite range.
The separate parity conjecture (1) is checked only against the published b-file
and is not a proposed settlement here.

`check.out` is the successful standard output; its final line is `ALL_CHECKS_OK`.
