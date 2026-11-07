# Rank-three webs in Gyotoku Example 15

Run `python3 docs/reports/gyotoku2607e15-rank-three-webs/check.py` with Python 3 and SymPy 1.14.0. The stored output ends with `ALL_CHECKS_OK`. The program uses exact symbolic expressions and integer arithmetic throughout.

The fidelity checks reproduce equation (18)'s Jacobian, Remark 14's three Gaussian relations, and footnote 4's inverse and beta–gamma density transformation in Yoshihiro Gyotoku, *The Independence-Preserving Property and Planar Web Geometry*, arXiv:2607.20646v1.

The primary proposed counterexamples are `F_lambda(x,y)=(x+y,lambda*x-y)` on the real plane, for positive lambda. The checks verify the inverse, nonzero Jacobian, three relations, the differential equation controlling all analytic relations, and the three-parameter solution. They also verify the chain-rule cancellation for the coordinatewise equivalence invariant `R=(u_x*v_y)/(u_y*v_x)`.

For each of the four families in Koudou–Vallois, *Independence properties of the Matsumoto–Yor type*, Bernoulli 18 (2012), Theorem 2.2, arXiv:1203.0381v1, the program computes two distinct exact values of R at points in the original positive quadrant. In the fourth family their difference is strictly negative for every positive delta by polynomial coefficient signs. Thus no listed map has the constant invariant of a Gaussian witness.

The differential-equation integration and the vector-space quotient dimension are justified in the preregistration's mathematical argument; printing their symbolic ingredients is not a proof of the dimension theorem. These checks are neither a Lean proof nor a complete replacement classification. The Gaussian and beta–gamma independence properties are classical; the proposed new comparison concerns only the conjectured exhaustive target list.
