# Reciprocal-sum square failures

`check.py` checks the proposed finite residue formula answering Section 7,
Question (1) of Matijasevick, Rodríguez and Salazar,
*Perfect squares in reciprocal-sum sequences and primes that are inert in
quadratic fields*, [arXiv:2610.00765v1](https://arxiv.org/abs/2610.00765v1).
The program uses only Python 3 standard-library exact integers and `Fraction`.
Run `python3 docs/reports/reciprocal-sum-square-failures/check.py` from the
repository root. Its deterministic output is stored in `output.txt`.

Fidelity checks reproduce Example 2.4 (h=3, x0=1) and the h=x0=10 example in
the proof of Theorem B. Direct rational arithmetic checks every h=1,...,80,
x0=1,...,160 and n=2,...,7 against the proposed failure indices and the
source's reduced-coordinate identity. Forward finite-field walks are checked
against independently constructed quadratic-preimage sets.

For every h=1,...,40 and h=50,65, every residue class modulo the proposed M,
and its shift by M, is checked against the complete eventually periodic
homogeneous-coordinate orbit modulo h. Cycles containing an obstruction are
rejected. The program compares the weighted product formula with all observed
classes and checks the quotient/remainder counting formula at seven cutoffs.

For h=65 it obtains M=4225 and counts (3507,714,4). Exactly the classes
1365,2015,2210,2860 have two failures; 520,1170,3055,3705 have coincident prime
obstructions and only one failure. At x0=1365 the first two failing readouts
are 484/13 and 214369/5.

These finite checks support a preregistration; they prove neither the
unbounded settlement nor the completeness of the literature screen. The
planned positive result is a general arithmetic and CRT argument, with these
computations serving as tests only. No Lean proof is part of this report.

The check program and README are independently authored for this report.
The source is cited for its definitions and fidelity examples; no source
TeX or third-party program is redistributed here.
