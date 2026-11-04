# Rectangular Mondrian defect eight

The source is Cooper O’Kuhn and Freeman Slaughter, *A Bounded-Defect Construction
for the Mondrian Art Problem*, arXiv:2610.19193v1, Example 1, Theorem 1 and
Conjecture 2. The check reproduces the published 315 by 315 example and the first
eight square-family members, and checks a proposed seven-tile rectangular family
on canvases of dimensions 15q by 42q with defect eight.

Run `python3 check.py`. The script uses only the Python 3 standard library and
exact integers and rational numbers. Its output is in `output.txt`.

The norm invariant and area identities are checked as exact polynomial
coefficient identities. Geometry is checked with pairwise interval separation,
an independent grid decomposition including boundaries, area sums and unordered
dimension pairs. Thirty-two recurrence members are checked, as are all norm-one
solutions with 1 ≤ q ≤ 1000. The affine map from the published norm-four family
is checked directly in every recurrence instance.

Finite computations do not establish the infinite cofinal statement. The
preregistration's proposed proof additionally uses induction for norm preservation,
positivity and unbounded growth. No Lean proof is included here, and no optimality
of eight is asserted. The program is an independent implementation of the cited
formulas and proposed construction.
