# Weighted convex-potential gap: fixed counterexample certificate

For Ahrami, arXiv:2605.24560v1, Open Problem 3.3, take

\[
\varepsilon=1/100,\qquad w(x)=1+\varepsilon(x-\pi/2),\qquad V(x)=\varepsilon x,
\quad 0\le x\le\pi.
\]

Both functions are affine. Thus the potential is convex, the weight is concave,
\(0\le V<1/25<2\), and \(49/50<w<51/50\).
The source parameters can be \(M=2,N_<=1/2,N_>=3/2\).

Run `python3 check.py` with Python 3 and SymPy. The recorded environment is
Python 3.14.8 and SymPy 1.14.0. Only integers, rational intervals and exact
symbolic algebra enter the evidence; no floating point calculations occur.
The expected final line is `ALL_CHECKS_OK`. Output is in `output.txt`.

The fidelity check reproduces the differential identity in the source's
Lemma 2.2. With \(q=\lambda w-V\), differentiating

\[
G(u'^2+qu^2)-G'uu'+G''u^2/2
\]

along \(u''=-qu\) gives \((2G'q+Gq'+G'''/2)u^2\).
At the Dirichlet endpoints the middle flux vanishes, yielding the printed
integral identity. Normalized sine eigenfunctions at \(V=0,w=1\), levels
1 and 2, give an additional exact check. The search seat's perturbation
coefficient is checked symbolically, but does not supply the fixed certificate.

## Certified shooting calculation

Put \(t=x/\pi\), and normalize the shooting solution by \(u(0)=0,u'(0)=1\)
where the prime now differentiates in \(t\). Its equation is

\[
u''=(A+Bt)u,\qquad
A=-\pi^2\lambda(1-\varepsilon\pi/2),\qquad
B=\varepsilon\pi^3(\delta-\lambda),
\]

with \(\delta=0\) for the zero potential and \(\delta=1\) for \(V=\varepsilon x\).
Writing \(u(t)=\sum c_kt^k\), the exact recurrence is

\[
c_0=0,\quad c_1=1,\quad
c_{k+2}=\frac{Ac_k+Bc_{k-1}}{(k+1)(k+2)},\quad c_{-1}=0.
\]

Machin's identity \(\pi=16\arctan(1/5)-4\arctan(1/239)\), with alternating
series remainders, encloses pi. Its angle identity is also checked using
integer real and imaginary parts of \((5+i)^4(239-i)\). Every interval
operation rounds outwards to the grid \(10^{-60}\mathbb Z\).

For every lambda tested, the program verifies \(|A|+2|B|<64\).
Replace coefficients by absolute values and evaluate at radius 2. The
Volterra terms are bounded by
\(2\,256^j/(2j+1)!\), so

\[
\sum_k |c_k|2^k\le 2\exp(16)<2\cdot3^{16}.
\]

Consequently truncation at degree 160 has error at \(t=1\) at most
\(R=2\cdot3^{16}/2^{161}<10^{-39}\). This estimate controls the full
shooting solution, rather than an unbounded perturbation remainder.
Adding \([-R,R]\) to the interval polynomial encloses the exact determinant
\(D(\lambda)=u(1)\). Bisection uses only strict certified signs, followed by
outward rounding and renewed sign checks at the delivered endpoints.

In units of \(10^{-12}\), the root brackets are:

| Potential | First eigenvalue | Second eigenvalue |
|---|---:|---:|
| 0 | [999989311923, 999989311927] | [4000050992847, 4000050992851] |
| x/100 | [1015707960628, 1015707960632] | [4015736948983, 4015736948987] |

## Continuum and ordering bridge

The analytic justification uses regular weighted Dirichlet Sturm–Liouville
theory, not a Galerkin truncation. For positive bounded weights and bounded
nonnegative potentials, the min-max of
\(\int(u'^2+Vu^2)/\int wu^2\) on \(H_0^1(0,\pi)\) gives the simple ordered
Dirichlet spectrum. For these continuous affine coefficients, a min-max
eigenfunction solves the classical ODE; uniqueness of the initial value
problem identifies its eigenvalues with the shooting zeros. Conversely a
shooting zero yields a nonzero eigenfunction. These continuum facts have
not been proved in Lean in this preregistration experiment.

For completeness, their mathematical construction is through the compact
positive self-adjoint operator \(T\) on \(H_0^1\), specified by
\(a(Tu,v)=\int wuv\), where \(a(u,v)=\int(u'v'+Vuv)\).
Poincare's inequality makes \(a\) coercive; Cauchy–Schwarz gives a uniform
one-half Holder bound for the unit energy ball, hence compact inclusion
into continuous functions by Arzela–Ascoli. Riesz representation and the
compact spectral theorem give \(T\); reciprocals of its positive eigenvalues
are the required min-max values. Weak eigenfunctions satisfy the ODE by
integration twice. The Wronskian proves simplicity. C1 Dirichlet trial
functions form a dense core, preserving each finite-dimensional min-max.
For the unweighted zero-potential operator the sine basis gives values
\(n^2\). Form comparison therefore gives, for either tested potential,

\[
\frac{50n^2}{51}\le\lambda_n\le\frac{50n^2+2}{49}.
\]

In particular, \(\lambda_1\le52/49\), \(\lambda_2\ge200/51\), and
\(\lambda_3\ge150/17\). The two sign-change brackets lie respectively
below \(200/51\) and between \(52/49\) and \(150/17\). Thus they enclose
exactly the first and second eigenvalues. This is why the four computed
root intervals certify spectral bounds.

Subtracting these enclosures gives the rigorous interval

\[
\Gamma[V,w]-\Gamma[0,w]\in
[-32692577,-32692561]/10^{12}<-3/100000<0.
\]

The planned Lean work is the continuum/min-max/shooting bridge and the
soundness of the rational interval and Taylor-tail certificate. This report
is exact external mathematical evidence, not a Lean proof of the claim or
its negation. It settles the displayed comparison with zero potential;
it makes no separate claim about minimization over all ordinary constants.
