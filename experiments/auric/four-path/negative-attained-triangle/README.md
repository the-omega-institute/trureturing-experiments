# Original negative attained-triangle certificate

The mathematical claim is an ordinary exact-arithmetic result. It has no Lean/kernel verification. Independent mathematical review and main-repository exposition are separate deliveries. The experiment inherits the repository's Apache-2.0 license; it embeds no third-party source code.

## Claim and original-model correspondence

Keep the original Dens2 x U(4) problem, ordinary transpose, the minimum of the two original sign distances, and a single common feedback scalar. Let

\[
h=1/\sqrt2,\quad s=\sin(\pi/8),\quad c=\cos(\pi/8),\quad
 e=2(h+s),\quad d=2(h+c),\quad H(\xi)=\xi/25+\xi^2/24.
\]

For every original pure effect and whole unitary with negative active sign and

\[
0<\rho(E,U)=\xi\le61/25,
\]

one has

\[
g(E,U)>H(\xi)>V(\xi).
\]

Here V is the attained positive geometric selector of published I37.5, not a global minimum. The negative population is q=(T^dagger E^T T)11, p=1-q. This statement includes every feasible population, either effect coherence phase, all attained-triangle boundaries, and x=1. It is stronger than the requested determinant-negative comparison. It implies no negative-active pure equality at or below the published selected zero, since IV49 gives xiStar<61/25. It asserts neither the complete two-sign firstloss/equality theorem nor completion of AURIC. No mixed-input pointwise inequality is claimed.

For a phase feedback the negative row is (a,b,ell,m)=(s,h,c,h). Put

\[
x=1-\xi/4+q/2,\quad L=1+x,\quad y=\sqrt{1-x^2},\quad
\mathcal T_x=\{(X,Y)\in[0,1]^2:\sqrt X+\sqrt Y\ge L\}.
\]

For a common point of that triangle set

\[
Z=\frac{L^2-X-Y}{2\sqrt{XY}},\quad
F=\sqrt{(1-X)(1-Y)}\sqrt{1-Z^2},
\]
\[
R^2=s^2X+h^2Y+sh(L^2-X-Y),\quad
C^2=c^2(1-X)+h^2(1-Y)+2chF,
\quad\Gamma=\sqrt p R+\sqrt q C.
\]

The exact triangle claim is e^2 p+d^2 q-4 Gamma^2>H(xi). Published II46.2 reconstructs R and C in one whole phase feedback; arbitrary effect phases have recovery modulus at most 2 Gamma. The arbitrary-unitary consumer is the finite same-effect, same-complex-active-trace transport of II42.3 and its nonnegative variance identity, used exactly as in II47.3. Because 61/25<13/4, every replacement has this same true minimum distance and unique negative active sign. A finite average of strict bounds at the same xi, plus nonnegative variance, remains strict. No inverse normal-orbit theorem is used.

## Upper bound for the same published V

Only the comparator lemma of the permitted sealed III50.4 is reused; its negative-canonical estimates and its pending canonical research are not premises. Its algebra and exact endpoint signs are independently checked in verify_fraction.py.

For 0<xi<=61/25 choose the unique zeta in (0,2/3) with

\[
\xi=\frac{4\zeta(4+3\zeta)}{2+3\zeta}.
\]

The map is increasing, with derivative 4(9 zeta^2+12 zeta+8)/(2+3 zeta)^2, and takes 2/3 to 4. The positive geometric point

\[
\bar x=1-\zeta,\quad \bar q=\frac{4\zeta}{2+3\zeta},\quad
\bar p=\frac{2-\zeta}{2+3\zeta},\quad \bar y=\sqrt{\zeta(2-\zeta)}
\]

has the same level xi and sqrt(qbar)*ybar=2 zeta sqrt(pbar). Its actual positive canonical value is

\[
G_+=\xi\frac{D+2hc\zeta+\zeta^2/2}{4+3\zeta},\qquad
D=4(2hs+s^2-hc).
\]

I37.5 gives V(xi)<=G_+. Direct subtraction gives

\[
\frac{H(\xi)-G_+}{\xi}
=\frac{P_H(\zeta)}{(4+3\zeta)(2+3\zeta)},
\]
\[
P_H(z)=(8/25-2D)+(18/25+8/3-3D-4hc)z+(84/25-6hc)z^2.
\]

Exact rational sqrt brackets verify that its quadratic coefficient is negative and both P_H(0) and P_H(2/3) are positive. Concavity then gives H(xi)>G_+>=V(xi). The positive comparison population qbar is distinct from the negative q; no independently optimized negative recoveries are combined.

## Small-x block

Feasibility and xi<=61/25 give x>=39/100 and q<=min(1,2x-39/50). On

\[
39/100\le x\le407/1000,
\]

put Qs=17/500 and define fixed positive constants

\[
R_s=s(407/1000)+h,\qquad
C_s=\sqrt{(c^2+h^2)(2-(139/100)^2/2)},
\]
\[
A_s=e^2-4R_s^2,\qquad B_s=d^2-4C_s^2.
\]

II46.21 gives R<=s x+h<=Rs. Since F<=sqrt((1-X)(1-Y)), Cauchy--Schwarz and X+Y>=L^2/2 give C<=Cs. Hence

\[
e^2p+d^2q-4\Gamma^2\ge \ell(q)
:=A_s+(B_s-A_s)q-8R_sC_s\sqrt{q(1-q)}.
\]

The rational verifier checks Rs>4/5, Cs>1, Bs-As<4 and ell(Qs)>2/5. For 0<q<=Qs<1/25,

\[
\ell'(q)=B_s-A_s-4R_sC_s\frac{1-2q}{\sqrt{q(1-q)}}<4-72/5<0.
\]

Continuity at zero gives ell(q)>=ell(Qs)>2/5, whereas H(xi)<=H(61/25)=1037/3000<2/5. No determinant or positive-population assumption is used.

## Large-x block: population chord and a joint matrix

Let t=1-x, with 0<t<=593/1000. The full triangle is parametrized by

\[
r=\sqrt X=1-ta,\quad z=\sqrt Y=1-tv,\quad v=(1-a)b,
\qquad 0\le a,b\le1.
\]

At a=1 one necessarily has z=1, so the redundant b coordinate includes that boundary. All other triangle points have these parameters by a=(1-r)/t and b=(1-z)/(t(1-a)). Define

\[
Q_p(t)=\min(1,61/50-2t),\quad h_0=4/25+2t/3,\quad H_0=t h_0,
\]
\[
\delta=2/25+2t/3+Q_p/6,
\quad A_0=e^2-H_0,\quad B_0=d^2-H_0-\delta.
\]

The two polynomial branches of delta are

\[
\delta=37/150+2t/3\ (0\le t\le11/100),\qquad
\delta=17/60+t/3\ (11/100\le t\le593/1000).
\]

For the actual feasible q in [0,Qp], direct expansion gives

\[
H_0+\delta q-H(4t+2q)=q(Q_p-q)/6\ge0.
\]

Put

\[
M=\begin{pmatrix}A_0-4R^2&-4RC\\-4RC&B_0-4C^2\end{pmatrix}.
\]

Then the original same-point lower objective minus H is at least
(sqrt(p),sqrt(q)) M (sqrt(p),sqrt(q))^T. The vector never vanishes, including q=0,1. Thus M positive definite is sufficient.

The exact quantities used by both verifiers are

\[
\mathfrak r=4sh-2s(h-s)a+2h(h-s)v
 -t[sh-s(h-s)a^2+h(h-s)v^2],
\]
\[
W=4\mathfrak r-h_0,\qquad
C_b=c^2a(2-ta)+h^2v(2-tv),
\]
\[
N=(1-a)^2ab(1-b)(2-ta)(2-tv)
 [4-t(1+a+v)][(2-t)^2-t^2(a-v)^2],
\]
\[
P=WB_0-4A_0C_b,\qquad
Q=r^2z^2P^2-16A_0^2c^2h^2tN.
\]

All factors of N are nonnegative on these boxes and r,z>=407/1000. The original joint formulas give exactly

\[
R^2=(h+s)^2-t\mathfrak r,\qquad
F=\frac{t^{3/2}\sqrt N}{2rz},\qquad
C^2=t\left(C_b+\frac{ch\sqrt{tN}}{rz}\right).
\]

Consequently

\[
M_{11}=tW,\quad
\frac{\det M}{t}=P-\frac{4A_0ch\sqrt{tN}}{rz},
\]
\[
Q=r^2z^2\left[P^2-
 \left(\frac{4A_0ch\sqrt{tN}}{rz}\right)^2\right].
\]

A0>4 follows from the exact constants and monotonicity of H0. The certificate proves W,B0,P,Q>0 on the full closed polynomial domain. Thus P exceeds the nonnegative radical, det M>0 and M11>0 for t>0. Sylvester's 2x2 criterion gives the required strict comparison. This is a common-implementation bound on R,C; their separate maxima are never substituted.

## Finite exact certificate and reproduction

certificate.json contains two complete binary partition trees with roots

- [0,11/100] x [0,1] x [0,1];
- [11/100,593/1000] x [0,1] x [0,1].

Coordinates are (t,a,b). Each internal node [j,left,right] bisects coordinate j at its rational midpoint. A leaf is 1. The 2385 leaves cover both closed roots, with shared boundaries included. They have 4768 total nodes and maximum depth 20. No sampled acceptance or discarded box is used.

verify_fraction.py is a standard-library-only checker. It constructs rational intervals for h,s,c using integer isqrt with mesh 10^-40. For each leaf it evaluates W,B0,P,Q by exact rational interval arithmetic. When needed, it intersects the natural interval with the multivariate mean-value interval f(m)+sum_j (box_j-m_j)*partial_j f(box), whose derivatives are computed by exact forward polynomial differentiation. Every test is a strict rational comparison. It also validates partition-node arity and coordinates. Its verified simple bounds on every leaf are:

| Polynomial | Strict rational lower bound |
| --- | --- |
| W | 12/5 |
| B0 | 49/5 |
| P | 1/10 |
| Q | 3/500 |

These are certificate bounds, not optimal extrema. verification.json retains the exact achieved rational minima, the small-x calculation, the comparator polynomial signs, and the known J-only method-failure enclosure. certify.py generates the finite tree using Arb at 128 bits; floating-point conversion chooses subdivision order only. The retained rational checker independently checks the tree and does not import Arb or the generator.

From this experiments worktree root:

```sh
python3 -B experiments/auric/four-path/negative-attained-triangle/verify_fraction.py experiments/auric/four-path/negative-attained-triangle/certificate.json experiments/auric/four-path/negative-attained-triangle/verification.json
python3 -B experiments/auric/four-path/negative-attained-triangle/certify.py --verify experiments/auric/four-path/negative-attained-triangle/certificate.json
```

Both commands exit 0 on the retained certificate. Python is 3.12.13; the second command requires python-flint 0.9.0. To regenerate the tree, use certify.py --seconds 300 --output PATH; a time-limited incomplete run never writes a complete certificate.

## Endpoint and original-unitary consumer

At x=1, II47.11 gives U=lambda Uminus, xi=2q>0 and g=d^2 q directly. Since q<=1, H(2q)<=37q/150<d^2 q. No coefficient containing 1/y is evaluated. The small-x and large-x blocks cover every other feasible x and include all pair, population, and distance endpoints through xi=61/25.

For an arbitrary whole U use published II42.3: every finite phase replacement has the same E, same complex negative active trace, same q,x and same true distance xi; g is their weighted mean plus nonnegative recovery variance. The strict phase bound therefore gives strict original g>H>V. Published III48.5 may then be used to exclude negative-active global minimizers on this prefix, but no pointwise purification of arbitrary mixed effects is inferred.

## Exact attained-gain corollary

For 0<x<1 define the original negative canonical amplitude and attained joint gain

\[
\Gamma_0=\sqrt p(sx+h)+\sqrt q\,c\sqrt{1-x^2},\quad
G_{\rm can}=e^2p+d^2q-4\Gamma_0^2,\quad
I_-=\max_{\mathcal T_x}(\Gamma-\Gamma_0).
\]

The maximum exists by II47.1: the same triangle is compact, XY>=x^2>0 and Gamma is continuous. Apply the proved all-point strict inequality at an actual maximizing triangle point. It gives exactly

\[
G_{\rm can}-V(\xi)-8\Gamma_0I_--4I_-^2
>H(\xi)-V(\xi)>0.
\]

The strict maximum-point argument includes every determinant-negative slice and does not replace Iminus by the failed Jminus envelope. No uniform positive margin as xi tends to zero is assumed.

## Sources and remaining scope

The fixed mathematical supplier is main commit 6b6539921bf3fde0f84f85fdda6c284b5891ae98: source23, I37.5, II42/46/47, III48.5, IV49. source_manifest.json seals their Git blobs and bytes. The permitted sealed comparator source is 4174415b64cd87a3770a85c30ec83cba6312cf7d:Continuation III, SHA256 2e595777683f204a62fdf5076591f0ee334ac069b279a84fbaeec4b256cc02e9, 44051 bytes. Its Lemma50.4 is the sole sealed mathematical input and its relevant algebra/signs are reproduced above. The negative-canonical gap, its independent in-flight task, and the positive-prefix candidate are not premises.

The mature ingredients are Cauchy--Schwarz, the quadratic chord identity, Sylvester's criterion, and interval inclusion/mean-value bounds. Vourdas arXiv:2203.09468v1, II.1(18), states a finite convex decomposition of a doubly stochastic matrix into permutations; published II42 consumes it only forward from an actual spectral decomposition. Jasper--Loreaux--Weiss arXiv:1602.08435v1 Theorem1.3 concerns prescribed singular values, not a fixed normal eigenvalue orbit, and supplies no converse used here. Boyd--Vandenberghe's convex optimization text does not supply an automatic exact SDP for this multi-radical triangle; no such equivalence is assumed. This result is a repository-derived synthesis, with no literature-wide novelty claim.

The unchanged II47 J-only completion route fails at xi=59/25,q=1/5,x=51/100. The rational checker reproduces all six supplied enclosures and Glow.upper<-1/5. This validates a known method failure; it is not an original E,U counterexample or an additional new result. The actual joined comparison is proved by the matrix certificate instead.

Complete two-sign ordering, full firstloss/equality, alpha/beta identification and the complete AURIC theory remain outside this result. Independent three-seat review and authorized main-paper retention remain pending. No push, PR, merge, or goal-state change is part of this experiment.
