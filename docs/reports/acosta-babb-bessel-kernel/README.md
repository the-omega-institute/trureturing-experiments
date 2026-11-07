# Acosta Babb's weighted Bessel kernel

`python3 check.py` uses only integer and `fractions.Fraction` arithmetic for
mathematical evidence. Its deterministic output is in `check-output.txt`.
The source is Ryan L. Acosta Babb, *A partial converse to the Riemann–Lebesgue
lemma for Bessel–Fourier series of order zero*, arXiv:2410.17681v1,
Section 3 (kernel and integral identity) and Section 4, Conjecture 1.
The journal version is *On the rate of growth of Bessel–Fourier coefficients
for integrable functions*, JMAA 559 (2026), 130440,
DOI:10.1016/j.jmaa.2026.130440, Sections 3 and 5.

The target is the explicit weighted display

    K_M^alpha(x) = sum_(m=1)^M (1-m/(M+1))*j_m^alpha*J_0(j_m*x), M>0,
    K_0^alpha(x) = 1.

The preceding printed S_m line contains j_m^alpha under its k-sum; it does
not imply the displayed Cesàro equality. This program uses the explicit
weighted display, which the source also uses in its integration arguments.

The fidelity checks reproduce the coefficients of the source identities
`integral_0^z y*J_0(y)dy = z*J_1(z)` and `J_0'=-J_1` through 40 terms.
Termwise integration/differentiation is justified by convergence of the
entire power series. The coefficients are the standard
`(-1)^k*(z/2)^(2*k+nu)/(k!*(k+nu)!)`, for nu=0,1.

The counterexample certificate uses 40 terms and an alternating tail. For
k>=40 the absolute-term ratios are below one and eventually below 1/2;
the first omitted term bounds the remainder. Rational endpoint signs give
roots a in (2.4048255,2.4048256) and b in (5.5200781,5.5200782).
To identify them as j_1 and j_2, pair consecutive series terms to see J_0>0
on [0,2]. For u(z)=sqrt(z)*J_0(z), the ODE is
`u''+(1+1/(4*z^2))*u=0`. On [2,6], comparison with the constant 17/16 gives
at least `4*pi/sqrt(17)>3` between distinct zeros. There is no root between
2 and a; an additional root between a and b would require total distance
larger than 6. These analytic arguments are not executed by Python.

The integral representation bounds |J_1'|<=1, so endpoint series bounds
transfer to the unknown roots and certify
`0.5191<J_1(a)<0.5192`, `0.3402<-J_1(b)<0.3404`.
Exact integer-power comparisons give `a^(7/5)<3.42`, `b^(7/5)>10.93`.
For every alpha in [2/5,1/2], factoring out `a^(alpha-2/5)` yields
`(K_2^alpha)'(1)>27843/500000>1/20`.
The integral representation gives |J_0''|<=1/2. Hence
`|(K_2^alpha)''|<895734319/60000000<15`.
Taylor's theorem at 1, where the kernel vanishes, then gives

    K_2^alpha(1-h) <= h*(-1/20+(15/2)*h) < 0,
    2/5 <= alpha <= 1/2, 0 < h <= 1/1000.

A separate direct interval evaluation certifies J_0(a*x), J_0(b*x) and the
kernel at alpha=1/2, x=999/1000, including interval uncertainty in a and b.
Square roots are enclosed by integer binary search, never floating point.

For the two-mode sign classification, put f(x)=J_0(a*x), g(x)=J_0(b*x).
On (a/b,1), f>0 and g<0, and the Wronskian identity implies that
`-g/f` increases strictly from zero to `b*|J_1(b)|/(a*J_1(a))`.
Thus the critical exponent is
`alpha_star=log(2*J_1(a)/|J_1(b)|)/log(b/a)-1`.
Exact integer-power comparisons enclose it between 17/50 and 7/20.
The Wronskian, zero-ordering and Taylor arguments are conventional analytic
proof obligations; this experiment does not claim Lean verification.
No classification for M>2 is asserted.
