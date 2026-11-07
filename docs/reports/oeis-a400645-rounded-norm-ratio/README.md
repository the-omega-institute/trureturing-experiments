# OEIS A400645: exact checks for the rounded norm ratio

Source: Bence Bernáth, [OEIS A400645](https://oeis.org/A400645), revision 21,
author line Oct 02 2026, revision timestamp 2026-10-11T22:07:30-04:00.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The recorded run used
Python 3.14.8. All computations use integers, rational numbers or exact
symbolic algebra; no floating point is used. Radical inequalities use rational
brackets certified by integer squares. The output is in `output.txt`.

The program reproduces all 51 published terms. Through index 4096 it checks
the signed radicands, independent nearest-integer midpoint brackets,
monotonicity, the Fibonacci upper bound, and the squared rounding-error bound.
It checks 263965 relative-growth inequalities, the quartic and octic root
identities, and the left-eigenvector identity used by the proposed proof.
The exact seed is a(26)=324; the proposed infinite-tail comparison has the
strictly positive rational margin 73872/17. Finite consecutive ratios at
indices 128 through 4095 are within 1/10^10 of the algebraic target, using
an exact isolating interval for that target.

These are finite checks and exact algebraic identities. They do not certify
the conjectured infinite limit or establish absence of prior literature.
