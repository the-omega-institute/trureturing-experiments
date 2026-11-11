# Perron-map differentiability checks

Source: Mykola Moroz, *Faithfulness and fractal (quasi-)equivalence principles
for Perron, Engel, and Pierce expansions*, arXiv:2510.03373v2, Sections 2, 3,
and the unnumbered question in Section 4.1 before Theorem 4.1.

Run `python3 check.py`. The output is in `output.txt`. Only integer arithmetic
and `fractions.Fraction` are used. The program is independently written from
the paper's mathematics; it incorporates no third-party program source.

The fidelity checks compare direct source-product sums with a separate
cylinder recurrence, reproduce the common positive/alternating diameter, and
evaluate constant-digit positive and alternating Lüroth series exactly.
The alternating summand includes its final factor `(q_k - 1)`.
The source's Section 3 adjacent-length ratio is printed reciprocally: its
diameter formula gives `(d+1)/(d-1)`, including the exact example `3` at `d=2`.

The settlement checks enumerate all depth 1 through 4 paths with each parent
parameter in 1 through 3 and each digit's excess over the minimum in 0 through
2. They check positive adjacency, alternating image order, the uniform
length ratio at most 3, the shrink bound, and signed secant bounds greater
than `1/8`. Independent coordinates test the geometric inequalities for all
sampled interval memberships, without asserting they specify the same code.
Five actual Lüroth points are additionally evaluated from eventually constant
codes at depths 1 through 16, with exact paired image values.

The general argument selects a point in the left adjacent cylinder whose
image is more than half a cylinder diameter away. The horizontal distance
is at most four such diameters. Alternating parity supplies both signs at
arbitrarily small scales. This is an unbounded mathematical argument;
the finite computations are checks and supply no formal proof of it.
