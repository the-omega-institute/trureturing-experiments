# Critical halo limit: exact preregistration checks

The target is the uncomputed critical limit following Theorem 6.6 of
Chon-Fai Kam, *Stellar rank under the contraction of SU(1,1) to the
Heisenberg–Weyl group*, arXiv:2610.00133v1, Appendix D.5.

Run with Python 3, using only its standard library:

    python3 check.py

`check.out` is the program's output. Every bound uses integers and
`fractions.Fraction`; no floating point is used as evidence.

The fidelity check reproduces the weight product in equation (5) and the
consecutive-term ratio in Appendix D.1 on 360 exact rational parameter
cases, including nonintegral Bargmann indices. It also checks both ratio
bounds used to control the head and infinite tail.

The settlement check encloses the positive critical root using the positive
series for the logarithm with a rational geometric remainder. The critical
identity gives lambda squared as `(1+3p)^3/(27p)`. Outward rational rounding
of the term recurrence and a geometric bound on the entire infinite tail
give certified enclosures for the critical sum, for integral and nonintegral
values of K. The cutoff always includes the outer hump; small terms never
trigger truncation. The proposed exact constant is enclosed by integer
square-root bounds and its prefactor identity is checked exactly.

The output brackets C between 3272073540231/10^12 and
3272073540232/10^12. It brackets S_1024 within 1/20 of C. Small values are
not monotone: S_2 is below C, and the tested sequence later overshoots it.
These finite checks support the proposed analytic evaluation; they do not
prove convergence for real K tending to infinity. The intended settlement
uses uniform Stirling asymptotics, endpoint estimates and Gaussian summation,
with finite checks serving only as helpers.
