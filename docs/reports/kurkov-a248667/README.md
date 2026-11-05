# Kurkov's A248667 composite-term factorization conjecture

`check.py` uses exact integers and `fractions.Fraction` to check the polynomial
family in [OEIS A248667](https://oeis.org/A248667), revision 53
(2026-03-12), and the proposed arithmetic route to Kurkov's conjecture.
Run from this directory with `python3 check.py`; `output.txt` is its output.
Python 3 and the standard library suffice. A failed assertion exits nonzero.

The script independently constructs each unreduced numerator by summing the
common denominator's polynomial tails. It reproduces the first six source
polynomials and all 60 displayed sequence terms. For indices 1 through 100 it
checks coefficient content, the constant and leading coefficients, and the
literal rational-expression identity at 0, 1 and 2. The leading coefficient is
n^(n-1), correcting the n^n typo in Taylor's comment.

For indices 1 through 3000 it checks the recurrence reset for 6846
prime-divisor pairs, membership by prime support, divisor closure, and internal
factorization of all 649 composite members. There are 1075 members in this
range. The script permits equal factors, including 9 = 3*3 and 49 = 7*7.
Membership above 100 is computed using Taylor's stated gcd criterion; it is
not an independent polynomial-content calculation in that range.

These finite checks support a proposed general argument. They neither prove
the unbounded conjecture nor constitute Lean kernel verification. No Lean
proof is part of this report. The program and output contain no copied source
files or third-party code; the embedded numerical examples are attributed to
the linked OEIS entry.
