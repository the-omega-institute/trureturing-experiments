# Fixed-field MRD column-limit check

`check.py` is a Python 3 exact-arithmetic check for the proposed limit in Gruica–Ravagnani–Sheekey–Zullo, Problem A / published Problem 2.6. It works over `F_2` using bit-vector Gaussian elimination and `fractions.Fraction`; no floating-point values enter any assertion.

It exactly counts ordered linear maps for `(q,n,d)=(2,2,2)` at `m=2,3`, independently enumerates all `[9 choose 6]_2 = 788035` subspaces for `(q,n,d,m)=(2,3,2,3)`, verifies the Gaussian exponents `G(2,1;2)=3` and `G(3,1;2)=7`, and brackets the proposed limits with a rational partial product and a union-bound tail estimate for `P(2)`. The captured stdout is `check.out`; the final line is `ALL_CHECKS_OK`.
