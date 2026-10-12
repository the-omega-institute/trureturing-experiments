# A square-variable theta derivative determinant obstruction

A computer-assisted Arb interval certificate gives, at both $T=0$ and
$T=1/10000$,

$$
\frac{\det[S^{(i+j)}(T)]_{i,j=0}^{4}}{S(T)^5}<-10^{25},
\qquad S(T)>0.
$$

The receiving kernel is the actual square-variable source
$S(T)=2K(\sqrt T)$. Thus the proposed full alternating derivative-Hankel
hierarchy fails at order five, where its required sign is positive.
The positive point is within the source's positive-half-line domain.
This is an obstruction to that stronger sufficient route, not to the
original prescribed-node scalar estimate, RH or Robin's inequality.
No Lean proof or literature originality is claimed.

## Source, variable and tested condition

Use the actual Csordas normalization

$$
K(u)=\sum_{m\ge1}\pi m^2(2\pi m^2e^{4u}-3)
                  e^{5u-\pi m^2e^{4u}},
$$

and Romik's original even analytic kernel

$$
\Phi(x)=2K(x/2),\qquad T=x^2/4,\qquad
S(T)=\Phi(2\sqrt T)=2K(\sqrt T).
$$

The source representation and full-kernel evenness are reused from
[Romik's orthogonal-expansion paper](https://doi.org/10.4064/aa200515-10-3),
equations (1.6)–(1.11). The square-variable kernel is also the one in
[Planat–Solé, arXiv:2608.19160v1](https://arxiv.org/abs/2608.19160v1),
equation (1.2), up to the positive factor two.
The preprint's Theorem 1.1 reports second-level concavity; its certificate
has not been rerun here and it does not assert the hierarchy tested below.

Write

$$
\Delta_k[S](T)=\det[S^{(i+j)}(T)]_{i,j=0}^{k-1}.
$$

The stronger proposed route requires

$$
(-1)^{k(k-1)/2}\Delta_k[S](T)\ge0
\qquad(k\ge1,\ T>0).
$$

For $k=5$ this sign factor is $+1$. The interval certificate directly
checks the opposite sign for this actual source. It does not reject the
order-four sign, second-level concavity or the scalar neighboring-row
node comparison. Original-variable derivative or Toeplitz counterexamples
are not used to transfer the negative sign to $T$.

For the actual positive exponential tilt $r_v(T)=A_ve^{T/v}S(T)$,
$A_v>0$, the standard binomial triangular derivative identity gives

$$
\Delta_k[r_v](T)=(A_ve^{T/v})^k\Delta_k[S](T).
$$

Consequently the same order-five hierarchy obstruction persists at each
fixed tilt $v>0$. This is a normalization consequence, not a prescribed
node inequality or a new derivative determinant theorem.

## Interval arithmetic and truncation

The program uses $K_{16}$, the sum of the first sixteen theta terms,
and derivative orders zero through eight. The full $K$ is even near
zero, whereas $K_{16}$ need not be even. At zero the exact full-source
conversion is

$$
S^{(j)}(0)=2j![u^{2j}]K(u).
$$

At a positive point the program computes the derivatives of the
finite evenization

$$
S_{16}(T)=K_{16}(\sqrt T)+K_{16}(-\sqrt T).
$$

Its Taylor expansion is holomorphic at $T=0$. The unsymmetrized finite sum $2K_{16}(\sqrt T)$ has odd powers of
$\sqrt T$ and is not used for derivatives at zero.
For each derivative the program adds the radius $2^{-180}$ derived below.
It then computes the determinant of the real interval matrix
$[S^{(i+j)}(T)/S(T)]_{i,j=0}^4$.
Its determinant is $\Delta_5[S](T)/S(T)^5$; the positivity of $S(T)$
and the strict upper bound $-10^{25}$ are interval checks.

## Analytic derivative remainder bound

On the complex circle $|u|=1/32$,

$$
\Re e^{4u}>1/2,\qquad |e^{4u}|<2,\qquad |e^{5u}|<2.
$$

For the first inequality use
$e^{-1/8}\cos(1/8)\ge(7/8)(127/128)>1/2$.
For the absolute bounds, $e^a\le(1-a)^{-1}$ for $0\le a<1$ suffices.
Using $3<\pi<4$ gives

$$
|K(u)-K_{16}(u)|
\le152\sum_{m\ge17}m^4e^{-3m^2/2}
\le304e^{-289}<2^{-280}.
$$

Here $m^4\le e^{m^2/2}$ for $m\ge17$, following from
$\log m\le m-1$ and $4(m-1)\le m^2/2$.
Writing $m=17+k$ yields $m^2\ge289+35k$ for integer $k\ge0$,
so the remaining Gaussian sum is less than twice its first term.
The series converges normally in a neighborhood of the circle, so
these Cauchy estimates apply to its analytic remainder.

Let $R=1/1024$. Evenizing the remainder defines a holomorphic function
$S-S_{16}$ on the $T$ disk $|T|\le R$, with

$$
|S(T)-S_{16}(T)|<2^{-279}.
$$

This follows by choosing either square root with $|u|=1/32$ on the
boundary and adding the two remainder bounds. The even power series
makes the definition independent of that choice and analytic at zero.
For $T_0\in\{0,1/10000\}$, Cauchy's estimate gives

$$
|S^{(j)}(T_0)-S_{16}^{(j)}(T_0)|
<\frac{j!2^{-279}}{(R-T_0)^j}
<2^{-180}\qquad(0\le j\le8).
$$

Indeed $R-1/10000>1/1200$ and the exact integer inequality
$8!\,1200^8<2^{99}$ bounds every order in that range.
The added radius dominates this analytic error; Arb series and matrix
arithmetic account for the finite rounding error separately.

## Reproduction and evidence scope

Requires Python 3.10+ and `python-flint==0.8.0`:

```sh
python3 -m venv /tmp/theta-hankel-venv
/tmp/theta-hankel-venv/bin/python -m pip install -r experiments/theta-square-hankel-obstruction/requirements.txt
/tmp/theta-hankel-venv/bin/python experiments/theta-square-hankel-obstruction/certificate.py --output /tmp/theta-square-hankel-result.json
```

The output path is explicit. `--bits` accepts at least 200 bits and
its default is 256. The script exits nonzero unless all point checks
pass. The retained result is the 256-bit output.

This certificate depends on Arb's interval series/determinant arithmetic
and the analytic Cauchy remainder bound. It supplies no Lean kernel
proof, no all-degree node comparison and no common-zero multiplicity
argument. Finite Jensen hyperbolicity already covered by complete
verified zero heights is reused rather than recertified here.
The all-integer Robin and RH goals remain open.
