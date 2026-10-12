# An actual-node obstruction to a cumulative Peano condition

An Arb interval computation encloses a degree-seven prescribed theta node
in

$$
0.03115806<v<0.03115808
$$

and checks a local sign that violates the proposed pointwise cumulative
Peano sufficient condition. This is a computer-assisted certificate with
analytic truncation bounds below; no Lean proof is delivered. It does not
refute the original neighboring-row estimate, RH, or Robin's inequality.
No claim of literature originality is made.

## Fixed source and receiving coordinate

Use the actual Csordas kernel

$$
K(u)=\sum_{m\ge1}\pi m^2(2\pi m^2e^{4u}-3)
                      e^{5u-\pi m^2e^{4u}},
$$

and the repository/Romik normalization

$$
\Phi(x)=2K(x/2),\qquad T=x^2/4,\qquad
S(T)=\Phi(2\sqrt T)=2K(\sqrt T).
$$

The receiving measure is fixed:

$$
d\mu(T)=\frac{2S(T)}{\xi(1/2)\sqrt T}\,dT.
$$

The source representation and its normalization are reused from
[Romik's orthogonal-expansion paper](https://doi.org/10.4064/aa200515-10-3),
equations (1.3), (1.6)–(1.11). The square-variable kernel is also the one
in [Planat–Solé, arXiv:2608.19160v1](https://arxiv.org/abs/2608.19160v1).
That preprint reports a finite-order concavity theorem; it does not assert
the cumulative condition tested here.

Put

$$
I_k=\int_0^\infty T^{k-1/2}S(T)\,dT
    =4\int_0^\infty u^{2k}K(u)\,du.
$$

For the monic Laguerre phase define

$$
\Pi_{n,v}(T)=\sum_{j=0}^n(-1)^{n+j}
 \frac{n!}{j!}\binom{n-1/2}{n-j}v^{n-j}T^j,
$$

$$
Q_n(v)=\int_0^\infty\Pi_{n,v}(T)T^{-1/2}S(T)\,dT.
$$

Multiplication by $2/\xi(1/2)>0$ preserves its zeros. With

$$
H_{n,v}(x)=\sum_{j=0}^n(-1)^j
 \frac{n^{\underline j}}{v^j}\frac{x^{2j}}{(2j)!},
$$

the coefficient identity is

$$
\Pi_{n,v}(T)=(-1)^n(1/2)_n v^n H_{n,v}(2\sqrt T).
$$

Here $(1/2)_n$ is a rising factorial. The identity uses
$(2j)!=4^j j!(1/2)_j$. Thus $Q_n(v)=0$ is the original prescribed-node
condition $v=n\lambda$ with `sourceQ n` vanishing at $\lambda>0$;
no new node scale or receiving measure is selected.

## Condition tested and the local obstruction

For a fixed prescribed node let

$$
\widetilde W(T)=T^{n-3/2}e^{-T/v},\qquad
\widetilde H(T)=\left(e^{T/v}S(T)\right)^{(n-1)},
$$

$$
D=\int_0^\infty\widetilde W\widetilde H,\qquad
Z=\int_0^\infty\widetilde W
  =v^{n-1/2}\Gamma(n-1/2).
$$

The proposed sufficient condition, when $D\ne0$, was

$$
\operatorname{sgn}(D)\mathcal P(t)\le0\quad(t\ge0),\qquad
\mathcal P(t)=\int_0^t(t-T)\widetilde W(T)
                    \left(\widetilde H(T)-D/Z\right)\,dT.
$$

Restoring the original positive constants in $g_v$ and $r_v$ multiplies
this signed expression by a positive factor and preserves its sign.
The condition is stronger than the scalar neighboring-row estimate;
its failure does not determine that estimate's sign.

For $n=7$, the certificate checks opposite signs of $Q_7$ at the rational
bracket endpoints and, throughout the whole bracket,

$$
D<-1/8,\qquad \widetilde H(0)-D/Z<-60{,}000{,}000,
$$

$$
D\left(\widetilde H(0)-D/Z\right)>9{,}000{,}000.
$$

Polynomial continuity gives at least one prescribed node in the bracket.
Uniqueness and a node's ordinal index are not required or claimed.
At any such node, continuity of $\widetilde H$ gives
$\widetilde H(T)-D/Z<0$ on some interval $[0,\varepsilon)$.
Because $\widetilde W(T)>0$ for $T>0$, one then has

$$
\operatorname{sgn}(D)\mathcal P(t)>0\qquad(0<t<\varepsilon).
$$

The endpoint values are calculated using the exact formulas

$$
D=\sum_{j=0}^6\binom6j v^{-(6-j)}(-1)^j
                  (11/2)^{\underline j}I_{6-j},
$$

$$
\widetilde H(0)=\sum_{j=0}^6\binom6j v^{-(6-j)}S^{(j)}(0),\qquad
S^{(j)}(0)=2j!\,[u^{2j}]K(u).
$$

The first formula follows by integration by parts in the rapidly decaying
actual theta series. At zero the required powers have positive exponents,
so the boundary terms vanish. The last formula uses the even analytic
kernel near zero. These formulas and the sign argument are ordinary
mathematical correspondences, not additional Lean declarations.

## Analytic bounds for omitted terms

The program uses $m\le16$, integrates $0\le u\le2$, and computes Taylor
coefficients through order twelve. It adds a radius $2^{-500}$ to each
$I_k$ for $0\le k\le7$, and $2^{-200}$ to each $S^{(j)}(0)$ for
$0\le j\le6$.

On $0\le u\le2$, using $3<\pi<4$,

$$
K_{>16}(u)\le32e^{18}\sum_{m\ge17}m^4e^{-3m^2}
 \le32e^{18}\sum_{m\ge17}e^{-2m^2}<2^{-536}.
$$

For $m\ge17$, $m^4\le e^{m^2}$; writing $m=17+k$ gives
$m^2\ge289+35k$, so the last sum is less than $2e^{-578}$.
Also $32e^{18}<2^{41}$ and $e>2$. The omitted contribution to each
moment on $[0,2]$ is at most $8\cdot2^{14}\cdot2^{-536}=2^{-519}$.

For $u\ge2$, $\log u\le u/2$ and summing the full theta series gives

$$
4u^{2k}K(u)\le256e\exp(16u-3e^{4u})\qquad(0\le k\le7).
$$

Indeed $m^4\le e^{m^2}$ for positive integers $m$, and
$\sum_{m\ge1}e^{-(3e^{4u}-1)m^2}<2e^{-(3e^{4u}-1)}$.
The tangent bound $e^{4u}\ge e^8+4e^8(u-2)$, together with $e^8>256$,
then gives an integrated tail less than $2^{-726}$. Combining both tails
is smaller than the added moment radius $2^{-500}$.

For Taylor coefficients, on the complex circle $|z|=1/32$,

$$
\Re e^{4z}>1/2,\qquad |e^{4z}|<2,\qquad |e^{5z}|<2.
$$

Hence the omitted analytic series is uniformly bounded there by

$$
|K_{>16}(z)|\le152\sum_{m\ge17}m^4e^{-3m^2/2}
 \le304e^{-289}<2^{-280}.
$$

Here $m^4\le e^{m^2/2}$ for $m\ge17$ and the sum of the resulting
Gaussian tail is less than twice its first term. Cauchy's coefficient
estimate gives

$$
|[u^{2j}]K_{>16}(u)|<2^{-280}32^{2j},\qquad
|S^{(j)}_{>16}(0)|<2^{-209}<2^{-200}\quad(0\le j\le6).
$$

The finite summands are entire, so the Arb complex-integral callback has
no unhandled branch cut. The change $T=u^2$ removes the moment endpoint
singularity. Finite Taylor coefficients are computed by Arb series
arithmetic, without differentiating a numerical remainder.

## Reproduction and evidence scope

Requires Python 3.10+ and `python-flint==0.8.0`:

```sh
python3 -m venv /tmp/theta-peano-venv
/tmp/theta-peano-venv/bin/python -m pip install -r experiments/theta-peano-node-obstruction/requirements.txt
/tmp/theta-peano-venv/bin/python experiments/theta-peano-node-obstruction/certificate.py --output /tmp/theta-peano-result.json
```

The output path is explicit. `--bits` accepts at least 180 bits and defaults
to 200. The script exits nonzero unless every interval sign check succeeds.
The retained `result.json` is the 200-bit output. The implementation was
run on CPython 3.13 with python-flint 0.8.0 on macOS; other platforms have
not been tested.

The certificate depends on Arb's interval arithmetic/integration and the
analytic tail arguments above. It is not a Lean kernel proof, a proof of
the all-degree neighboring-row estimate, or a common-zero multiplicity
argument. The proof obligations for RH and strict Robin at every integer
greater than 5040 remain open.
