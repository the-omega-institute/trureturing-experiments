# Simple random walk Hardy weight: fourth-order coefficient

The target is the unnumbered conjectural Remark after Proposition 10 in Hake, Keller and Pogorzelski, *Optimal Hardy Inequalities for Random Walks on Z²*, arXiv:2608.24213v1, Section 4. The operator is I-P with probability 1/(2d) at every signed coordinate step.

Run `python3 check.py` with Python 3.12 or later, SymPy 1.14 and mpmath 1.3. The exact part alone is available with `python3 check.py --exact-only`.

The exact part uses symbolic rational algebra. It reproduces the source Lemma 11 covariance identity for simple random walk in dimensions 3 through 8, symbolically in all coordinates. It verifies the reciprocal Fourier-symbol expansion, harmonic quartic, radial fourth derivatives, and the identity between the resulting Hardy correction and the source polynomial Q. It also checks 191 nonzero small lattice directions up to permutation symmetry, in dimensions 3 through 6. These checks do not prove the asymptotic estimates.

The numerical part is explicitly separate from the exact evidence. It evaluates the independent heat-kernel integral G(x) = integral from 0 to infinity of exp(-t) times the product of I_|x_j|(t/d). It uses 42 decimal digits of working precision, tanh-sinh quadrature, and a 40-term scaled-Bessel expansion only when its argument exceeds both 10000 and 100(|x_j|+1)². At each sampled point it checks relative Green harmonicity to 10^-25. It computes r^6 times the conjectured expansion residual in dimensions 3 and 4 along the axis, a face diagonal, and the full diagonal, with scales 4, 8 and 12. The numerical assertions only check finiteness and a sample bound of 100. They provide no certified numerical error or bound for unsampled lattice points.

`check.out` is the program's output. Its final line is `ALL_CHECKS_OK`. No Lean proof or negation of the conjecture is supplied by this report.
