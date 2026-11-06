# Alm–Maddux RA52 Boolean-group family

Exact finite checks for Problem 3 of Jeremy F. Alm and Roger D. Maddux,
*Finite representations for two small relation algebras*,
[arXiv:1712.00129v1](https://arxiv.org/abs/1712.00129v1),
submitted 30 November 2017; Algebra Universalis 79, article 87 (2018),
[DOI 10.1007/s00012-018-0570-4](https://doi.org/10.1007/s00012-018-0570-4).

Run with Python 3 (standard library only):

    python3 docs/reports/alm-maddux-ra52-boolean-family/check.py

The program reproduces the source's explicit dimension-ten representation
from its 63 listed subgroup bitstrings, checks all six atomic sumsets by
direct pair enumeration, and reproduces its Johnson-scheme point count and
strict rational union bound for n = 13. The transcribed numerical data are
credited to Alm and Maddux; no source program is copied.

The proposed family uses k triples and one final coordinate. H has final bit
zero and even parity in every triple. Set B = H minus zero, C = {wt > 2k},
and A = G minus (H union C). Integer Walsh–Hadamard XOR convolutions check
all six sumsets at k = 3, 4, 5, with a direct pair-enumeration cross-check at
k = 3. The source cycle list independently reconstructs the table. The
binomial recurrence and strict density bound are checked for 2 ≤ k ≤ 100.

All expected family identities pass for k = 3, 4, 5. At k = 2 the only
failure is B+C missing integer 63 (the six triple coordinates are one and
the final coordinate zero). This limits this construction; it does not
refute the source's dimension-seven Problem 2. The recorded output ends
with ALL_CHECKS_OK and the program exits zero. Finite tests do not establish
the universal theorem, and no Lean proof or kernel verification is claimed.

Coordinates in the proposed family use bit i as coordinate i, with final
bit 3k. Source bitstrings are parsed in their printed MSB-first order.
All evidential arithmetic is exact: Python integers and fractions.Fraction.
