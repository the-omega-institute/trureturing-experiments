# Hanna's partition-multinomial exponential integrality

The two targets are Paul D. Hanna's conjectures in [OEIS A183241](https://oeis.org/A183241)
(revision 12) and [OEIS A182963](https://oeis.org/A182963) (revision 8). Both assert
integrality of the coefficients of exp(sum b_n x^n/n), with b_n the partition-indexed
sum of squared or cubed multinomial coefficients, respectively. Their definitions
are supplied by A183240, A183235 and A036038. Each integer partition contributes
once; repeated equal parts introduce no additional factorial denominator.

Run `python3 check.py`. Python 3's standard library suffices; every calculation
uses integers or `fractions.Fraction`. The stored result is `output.txt`.

The source coefficients are reproduced through degree 15 for A183241 and degree
13 for A182963, together with both logarithmic inputs through degree 7. For each
power r from 1 through 5 and degree 0 through 20, independent partition sums and
weighted partition products agree, and the rational recurrence, direct formal
exponential and integral Euler product agree with nonnegative coefficients.
The divisor identity has nonnegative integral primitive-orbit counts throughout
that range. Direct word enumeration additionally checks cardinalities, every
least-period root/repetition bijection and actual primitive orbit sizes for
r = 1 through degree 6, r = 2 through degree 5, and r = 3 through degree 4.

These are finite fidelity and consistency results. They do not establish an
all-index theorem or provide a Lean proof. The proposed general mechanism is
simultaneous cyclic rotation of equal-content word tuples, with positive weakly
decreasing letter multiplicities retained under taking a repetition root.
