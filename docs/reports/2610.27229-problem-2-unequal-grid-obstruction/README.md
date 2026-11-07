# Cournot letter-H obstruction: exact checks

Source: Juan Luis García Guirao, *Connected ω-limit sets of Cournot maps: periodic orbits, realization and almost every point*, arXiv:2610.27229v1, Section 11, Problem 2. The source asks whether H = ({0,1} × [0,1]) ∪ ([0,1] × {1/2}) is a point ω-limit set of a continuous Cournot map F(x,y) = (f₂(y),f₁(x)).

Run `python3 check.py`. The program uses Python integers and `fractions.Fraction` exclusively. The companion `check-output.txt` ends with `ALL_CHECKS_OK`.

The fidelity checks reproduce the source's even-iterate identity for 729 continuous piecewise linear reaction-function pairs, nine starting points, and n = 0,…,6 (45,927 exact equalities). They also check the displayed map in the source's robust-lattice example: 84 exact contraction identities and its exact image of [1/2,1].

For every pair of continuous piecewise linear maps whose knots and knot values belong to the three-point or four-point uniform rational grid, the program constructs a point of H absent from F(H). These are 729 and 65,536 pairs respectively. Membership in F(H) is evaluated from the exact three segment images, using full interval ranges rather than sampled points. The candidate ordinate avoids the two values f₁(0), f₁(1), and its abscissa differs from f₂(1/2).

For 1 ≤ |P|, |Q| ≤ 5, the program exhausts 11,398 finite function assignments to check the two necessary surjection constraints for a strongly invariant rectangular grid. They are compatible exactly when the two cardinalities coincide; the diagonal counts are (m!)². This checks the proposed finite obstruction and does not supply the balanced-grid realization theorem.

These finite checks do not certify the theorem for arbitrary continuous maps, existence of dense orbits, or strong invariance of arbitrary ω-limit sets. A general argument using strong invariance is required to settle Problem 2. No Lean proof is included.
