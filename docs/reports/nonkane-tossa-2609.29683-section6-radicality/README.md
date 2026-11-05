# Nonkané–Tossa Definition 3.1: exact radicality checks

Source: arXiv:2609.29683v1, Section 6, fourth further-direction bullet.
This experiment tests the literal polynomial definition, not a representation
theoretic replacement.

Run `python3 check.py` with Singular 4.4.1 and its bundled `primdec.lib`.
The Python standard library enumerates all partitions, bijective tableaux and
supports. Singular constructs the source ideals by ideal sums and intersections,
computes `radical(I)`, and compares ideals by exact Gröbner containment in both
directions. It uses rational or algebraic-number coefficients and prime fields;
there is no floating point evidence. Output is in `output.txt`.

The range is r in {2,3}, n in {0,1,2,3,4}, every z in {0,...,n}, and every
partition of n-z: 104 cases, with 82 radical and 22 nonradical ideals. Every case
agrees with radicality iff r is nonzero in the coefficient field or the partition
has at most one row. The empty partition gives the unit ideal.

Characteristic-zero fields are Q(zeta_2)=Q and Q(zeta_3), the latter represented
by a with a^2+a+1=0. The inseparable cases use F_2 for r=2 and F_3 for r=3.
Those finite fields test the polynomial formulas; they do not satisfy the paper's
standing infiniteness hypothesis. The same explicit nilpotents persist over
an infinite extension of the prime field.

Fidelity checks reproduce Examples 4.1 and 4.3 at n=3, z=1, mu=(1,1):

    (x1*x2*(x1^r-x2^r), x1*x3*(x1^r-x3^r),
     x2*x3*(x2^r-x3^r), x1*x2*x3).

For r=p=2 and r=p=3 the exact ideal reductions verify that
w=x1*x2*(x1-x2) is outside this ideal and w^p is inside it. These are global
nilpotents, and the test does not infer nonradicality from point samples.

The finite range is evidence for a proposed general argument, not a proof of
the unrestricted characterization. No Lean proof is part of this report.
