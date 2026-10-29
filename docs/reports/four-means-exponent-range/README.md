# Four means: real exponent range

Exact checks for Romeo Meštrović and Miomir Andjić, *Two curious inequalities involving different means of two arguments*, arXiv:1804.00542v1, Remark 1.5 and Conjecture 1.6; Resonance 28 (2023), 1891–1897, DOI 10.1007/s12045-023-1726-x (journal Conjecture 1.5).

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The recorded stdout is `check.out`. Success requires exit 0 and final line `ALL_CHECKS_OK`. No floating-point calculation is used as evidence.

The source defines H=2xy/(x+y), G=sqrt(xy), A=(x+y)/2 and Q=sqrt((x²+y²)/2), for positive x,y. Its conjecture asserts A^n+G^n≤Q^n+H^n universally for n<0 or n≥1/2, and noncomparability for 0<n<1/2.

The fidelity check reproduces Theorem 1.1's product inequality and its exponent-two inequality using symbolic identities and 484 positive rational pairs. The source's separate equality clause needs an exception at exponent zero; the conjecture does not include that clause.

For u=G/A, put F_n(u)=(2-u²)^(n/2)+u^(2n)-1-u^n. The check verifies F_n(1)=F'_n(1)=0 and F''_n(1)=4n(n−1), together with the coefficient n(n−1)/2 at ε⁴ for x=1+ε,y=1−ε. These identities support a local obstruction throughout 0<n<1; continuity and Taylor's theorem remain part of the mathematical argument, rather than consequences of a finite grid.

An infinite exact family at n=1/2 uses u=t², x=1+sqrt(1−t⁴), y=1−sqrt(1−t⁴), with 1/2≤t<1. The polynomial identity

    (1+t−t²)^4 − (2−t⁴) = (t−1)²(t+1)(t⁵−3t⁴+4t²+3t−1)

has a strictly positive right side: the last polynomial on [1/2,1] has positive Bernstein coefficients, all checked exactly. Both compared fourth roots are positive. Thus the conjectured inequality fails on this entire family. A separate algebraic anchor is n=1/2, x=2+sqrt(3), y=2−sqrt(3).

The 45 rational exponents p/q with 0<p<q≤12 are each tested at actual u∈(0,1) giving both strict signs. Rational bisection bounds each rational power; every root enclosure verifies lo^q≤target≤hi^q. These finite checks corroborate the family argument and do not certify the universal real-exponent classification or its literature novelty. The global algebraic identities used in the classification outline are also checked.

Source: <https://arxiv.org/abs/1804.00542v1>. The source is quoted only by reference; the check program is independent code.
