# An actual-node screen of a stronger Laguerre correlation condition

This experiment encloses one proposed sufficient condition on the unchanged
Riemann theta source. At its exact degree-one Jensen node, the condition
fails. This is an ordinary numerical enclosure, without Lean certification.
It gives no negative response sign, first defect, RH counterexample, or
conclusion about a restricted family of high-degree nodes.
It gives no uniform strict Robin conclusion for every $H>5040$;
RH and that universal Robin claim remain unresolved.

The numerical result is

$$
c^2+r^2-1\in[-0.290064,-0.156722]<0.
$$

The result excludes this condition as a universal requirement for every
actual node. It does not exclude weaker joint conditions. Cauchy--Schwarz,
the Laguerre Rodrigues formula, and the original theta representation
remain existing mathematical tools; the experiment does not claim them
as new results or certify bibliographic originality.

## Source and exact node

Use the original Romik normalization

$$
\Phi(x)=\sum_{m\ge1}
\left(4\pi^2m^4e^{9x/2}-6\pi m^2e^{5x/2}\right)
e^{-\pi m^2e^{2x}},\qquad x\ge0,
$$

with even extension. Write
$M_j=\int_0^\infty x^j\Phi(x)\,dx$ and set

$$
\alpha=M_0/M_2,\qquad t=2\alpha.
$$

For the normalized law $\Phi(x)\,dx/\Xi(0)$ on the whole line,
$B_{1,t}(x)=1-tx^2/2$ and $J_1(t)=\mathbb E B_{1,t}$.
Thus this $t$ is its unique actual node, enclosed by

$$
t\in43.28068807484977499728988899251696
\ \pm\ 3.61\times10^{-33}.
$$

Let $y=\alpha x^2/2$, $a=3/2$,
$U(y)=e^y\Phi(\sqrt{2y/\alpha})$, and

$$
f(y)=U'(y),\quad
H(x)=\Phi(x)+\frac{\Phi'(x)}{\alpha x},\quad
R(y)=\frac{\int_0^y s^{1/2}e^{-s}\,ds}{y^{1/2}e^{-y}}.
$$

Then $f=e^yH$. Define $\nu$ by normalizing
$|f(y)|y^{1/2}e^{-y}dy$, equivalently $|H(x)|x^2dx$.
For $S=\operatorname{sign}f$, $Y=y$, the screened condition is

$$
c^2+r^2>1,\qquad
c=\operatorname{Corr}_\nu(S,Y),\quad
r=\operatorname{Corr}_\nu(Y,R(Y)).
$$

The exact centering and signed moment follow by integration by parts:

$$
\int_0^\infty H(x)x^2dx=M_2-M_0/\alpha=0,
$$

$$
\int_0^\infty YH(x)x^2dx
=\alpha M_4/2-3M_2/2.
$$

Consequently $\mathbb E_\nu S=0$ and $\operatorname{Var}_\nu S=1$.
The program encloses the six raw integrals with weights
$1,Y,Y^2,R,R^2,YR$ against the same density, and uses these identities
for the signed moment. Outward interval arithmetic conservatively encloses
the shared-input dependencies; it does not explicitly preserve correlations.

## Enclosures and reproduction

The program reuses `polys` from
[derivative_bandwidth.py](../../docs/reports/theta-mixed-matrix/derivative_bandwidth.py)
through AST extraction, without executing that module's experiment.
It checks the supplier's SHA-256 before use.
The source representation is Romik, arXiv
[1902.06330v3](https://arxiv.org/abs/1902.06330v3), equations (1.6)--(1.11).
The repository's exact finite test normalization is recorded in
[LC423--426](https://github.com/the-omega-institute/trureturing/blob/360dae7811de072f5d31801d99c50ffb4b6a96a6/Library/Analytic/romik2021orthogonal.md).

The retained theta terms are $m=1,\ldots,8$. All $m\ge9$ terms are
enclosed by decreasing-ratio geometric majorants, including the first
two derivatives. The source moments use Arb analytic quadrature on
$[0,3]$ and explicit omitted-source bounds. For the six absolute-density
moments, 131,072 interval rectangles cover $[1/4096,3]$ at 160 bits.

The tiny core uses $\Phi'(0)=0$ and the enclosed $\Phi''$ to bound
$|H|$; its raw mass allowance is below $8.107\times10^{-12}$.
For $R$, the program uses

$$
R(y)=\frac{\sqrt\pi\,e^y\operatorname{erf}(\sqrt y)}{2\sqrt y}-1
=y\int_0^1v^{1/2}e^{y(1-v)}dv,
$$

which is increasing and satisfies $R(y)\le(2/3)ye^y$.
The enclosed node has $20<\alpha<22$.
For every $x\ge3$, the full-source envelopes
$|\Phi(x)|\le64e^{9x/2-\pi e^{2x}}$ and
$|\Phi'(x)|\le800e^{13x/2-\pi e^{2x}}$,
together with $R(y)\le e^y$ for $y\ge1$, bound all six raw integrands by

$$
121\cdot864\,e^{-2e^{2x}}.
$$

Here $x^6\le e^{6x}$, $22x^2\le e^{2x}/2$, and
$(25/2)x\le e^{2x}/8$ on that range. Integrating the displayed
envelope gives a bound below $2.50\times10^{-349}<2^{-1000}$.
The source moments have the same spatial tail allowance.
No exterior cutoff is treated as the complete source.

Use Python with `python-flint==0.9.0`:

```sh
python experiments/jensen-laguerre-correlation-screen/n1_certificate.py
```

[results.json](results.json) stores the numerical enclosures, source hash,
runtime, and exact dyadic endpoints of the score interval. The program
fails unless the final interval certifies a strictly negative score.
