# Actual c27 prime debt and even high-frequency floor

The target is the complete Weil form on `(-a,a)`, `a=3log(3)/2`, with
all fourteen prime powers below 27. The endpoint shift at 27 is zero in
L2. This packet supplies a smaller prime debt for the existing even
Neumann–Gamma comparison. The retained matrix sign and its coupling to
the entire complement are not certified.

The weighted Schur formula is reused from Liu, *Certified Weil Positivity
Beyond the Unit Window*, Appendix B.4, pinned external commit
`b6cd2183c1e79c6c27a34267812a7b2d73ed1b59`. Its old support constants are
not reused. The earlier project application at c9 is in
[`Library/Weil/liu2026tailcompensation.md`](https://github.com/the-omega-institute/trureturing/blob/c1427c902273fc8a9269cd2e0e07fae450a9d60f/Library/Weil/liu2026tailcompensation.md).
The Neumann–Gamma high-complement comparison is already in the project's
`UNIFORM_SCHUR_CERTIFICATE_ACROSS_PRIME_THREE` section, PT7–PT10, and its
periodic Gamma/form-domain suppliers G1–G3. No new general Schur,
Neumann or spectral approximation theorem is supplied.

For the actual clipped prime translation sum `C_p`, set
`q(u)=1+(29/50)(u/a)^2`. Exact Fraction minima on every outward support
cell certify `||C_p|| <= 59/10`. Every quadratic margin exceeds `1/400`;
the report retains the stronger exact minimum. Rational logarithm
enclosures use a `10^12` grid from directed interval logarithms; square
root lower bounds use checked integer squares. The coincident support
switches from 3 and 9 use their exact normalized shifts `2/3` and `4/3`.
No mesh is substituted for an interval supremum.

The sum of the individual paired-shift norms gives the larger debt
`8.54672965014089...`. With that debt, the first positive bound
`A(2pi(N+1)/log27)-debt` occurs at N=16973. With the conservative joint
debt `59/10`, it occurs at N=1203. Here
`A(t)=Re digamma(1/4+it/2)-log pi`.
At N=1280 the latter floor lies in
`[0.06273617862538,0.06273617862539]`, strictly above `1/16`.
These are thresholds for the specified sufficient high-complement
estimate, not necessary matrix sizes or positivity results for the full
form. A negative value of this lower bound is not a negative Weil test.

On the existing paper-level even Neumann–Gamma and form-domain
identifications, the scalar packet therefore supplies
`Q(f)>=||f||^2/16` for every even f perpendicular to the cosine modes
0 through 1280. Even pole energy remains nonnegative. No low-block sign,
full Schur certificate, Lean proof or cofinal-support result is claimed.

## Reproduction

The existing arithmetic module is read-only input, fetched from the main
project's fixed commit. Its SHA256 is
`8bb067fc5499b0f2e1e48836e7a82237a15504109f82a856c72478d1096d69d0`.
Only `psi_iv` is called for the Gamma values; its old fixed-window run is
not executed. The arithmetic module is not redistributed here.

From this repository root, with the main project at `../trureturing`:

```sh
python3 -m venv build/weil-c27-venv
build/weil-c27-venv/bin/python -m pip install -r docs/reports/weil-c27/requirements.txt
mkdir -p build/weil-c27-inputs
git -C ../trureturing show 4804c4020d7d9f165a2e2683c42b01e3d47a8be4:research/weil_ground_mode/certify_prime3_refined.py > build/weil-c27-inputs/certify_prime3_refined.py
build/weil-c27-venv/bin/python docs/reports/weil-c27/c27_high_floor.py --arithmetic-source build/weil-c27-inputs/certify_prime3_refined.py --out build/weil-c27-high-floor.json
```

The checked run uses CPython 3.9.6, mpmath 1.3.0, NumPy 2.0.2,
SymPy 1.14.0 and 70 interval digits.
The pinned NumPy and SymPy dependencies support importing the existing
arithmetic module; its array/certificate routines are not run. The
machine result is `high-floor.json`; it is numerical input evidence,
not a proof-assistant or independent reproduction claim.

## Precision-parametric arithmetic input

`arb_arithmetic.py` applies the existing actual Fourier assembly to c27
using FLINT's `acb.digamma` and `acb.polygamma(1)`. These are existing
certified special-function suppliers, also used elsewhere in this
experiment repository; their algorithms are not reimplemented here.
The interface replaces the old fixed special-function remainder only.
It does not rerun the old spectral certificate or change its parameters.

For full length $d=\log c$, write $\omega_n=2\pi n/d$,
$z_n=1/4+i\omega_n/2$, $b_j=2j+1/2$, $H=\cosh(d/2)$ and
$w_r=\Lambda(r)/\sqrt r$. The scalar formulas are the existing
[arithmetic boundary symbol](https://github.com/the-omega-institute/trureturing/blob/212cac9d3e065e29e07f43b0ec5e48b2a066640c/D5/S3/Weil/ZetaBridge/WeilArithmeticCouplingJet.lean)
and [actual Fourier assembly](https://github.com/the-omega-institute/trureturing/blob/4804c4020d7d9f165a2e2683c42b01e3d47a8be4/research/weil_ground_mode/certify_prime3_refined.py),
with every prime power $r<c$ retained:

$$
\begin{aligned}
s_n={}&-\frac{2\omega_n(H-1)}{\omega_n^2+1/4}
-\frac12\Im\psi(z_n)
+\sum_{j\ge0}\frac{\omega_nc^{-b_j}}{b_j^2+\omega_n^2}
-\sum_{r<c}w_r\sin(\omega_n\log r),\\
a_n={}&\Re\psi(z_n)-\log\pi+
\frac{\Re\psi_1(z_n)}{2d}
-\frac2d\sum_{j\ge0}c^{-b_j}\Re(b_j-i\omega_n)^{-2}\\
&+\frac{4(H-1)}d\Re(1/2+i\omega_n)^{-2}
-2\sum_{r<c}w_r(1-\log r/d)\cos(\omega_n\log r).
\end{aligned}
$$

The endpoint shift at integer c is zero almost everywhere. The actual
unphased orthonormal even basis is $1/\sqrt d$ and
$\sqrt{2/d}\cos(2\pi nx/d)$ on $(-d/2,d/2)$. Its entries are

$$
E_{00}=a_0,\qquad E_{0n}=-(-1)^n\frac{\sqrt2s_n}{\pi n},\qquad
E_{nn}=a_n-\frac{s_n}{\pi n},\qquad
E_{mn}=(-1)^{m+n}\frac{2(ns_n-ms_m)}{\pi(m^2-n^2)}\quad(m,n>0,\ m\ne n).
$$

`Arithmetic.even_entry` consumes the supplied scalars in these formulas.
The original divided-difference formula without these parity factors
uses the phase-adjusted basis $(-1)^n\sqrt{2/d}\cos(2\pi nx/d)$.
The diagonal congruence transports that formula to the unphased basis
declared by this interface, including its zero column.

The two decaying Gamma corrections are truncated after R terms. Set
$q=c^{-2}$ and $b_R=2R+1/2$. Uniformly over all nonnegative frequencies,
the omitted absolute contributions are bounded by

$$
\epsilon_s=\frac{c^{-b_R}}{2b_R(1-q)},\qquad
\epsilon_a=\frac{2c^{-b_R}}{d b_R^2(1-q)}.
$$

Here $\omega/(b^2+\omega^2)\le1/(2b)$ and
$|(b-i\omega)^{-2}|\le1/b^2$; subsequent terms have geometric
decay q and nondecreasing denominators. The program adds these explicit
tail balls to the special-function and finite-arithmetic balls. It
chooses R from the requested output precision instead of retaining a
fixed asymptotic remainder. The exact zero symbol $s_0=0$ is preserved.

The default packet supplies all modes 0 through 1280 and three additional
coupling samples, 1281, 32768 and 32769. It uses 768 working bits and a
$2^{-640}$ output grid. Each scalar interval has width at most two grid
units; its midpoint therefore has error at most $2^{-640}$. This is a
scalar accuracy statement. Subsequent matrix operations, infinite
coupling sums and a positivity decision need their own error bounds.

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/arb_arithmetic.py --extra-mode 1281 --extra-mode 32768 --extra-mode 32769 --out build/c27-arithmetic.json.zlib
```

The output is zlib-compressed JSON. `rows` uses exact integer endpoint
pairs on the `grid_bits` dyadic grid; the two tail-bound fields use
`tail_grid_bits`. Only the scalar table is supplied: no matrix sign,
full infinite-tail certificate, Lean verification or RH conclusion is
asserted. FLINT's documented APIs are
[digamma/polygamma](https://python-flint.readthedocs.io/en/latest/acb.html)
and its [ball-arithmetic contract](https://flintlib.org/doc/using.html).

The saved `arithmetic-scalars.json.zlib` packet contains 1,284 rows. All
rows agree exactly on the dyadic output grid between the 768-bit and
1024-bit executions; the observed maximum interval width is one grid
unit. Eight c3 diagonal samples and seven nonzero c3 symbol samples are
contained in the older read-only scalar routines' intervals. Those calls
check interface normalization only; the older certificate was not run.
The retained packet SHA-256 is
`f993be1aa48824cbf22276ac2d80e834ca9e54e4519e2040f437ebb4f724990f`.

## Infinite moments of the same actual symbol

`arb_tail_moments.py` supplies directed intervals for

$$
\mathsf T_p=M^{p-1}\sum_{m>M}\frac{s_m}{m^p},\qquad
\mathsf U_p=M^{p-1}\sum_{m>M}\frac{s_m^2}{m^p},\qquad p\ge2.
$$

These are discrete moments of the symbol defined above, including the
pole, full Gamma contribution and all fourteen c27 prime powers. They
are different from the integrated Gamma moments used by the c9 source.
The normalization prevents powers of the exterior cutoff from obscuring
the accuracy needed by a coupling consumer. Divide a supplied moment by
$M^{p-1}$ to obtain its unnormalized counterpart.

The reused analytic input is Fredrik Johansson,
[*Arbitrary-precision computation of the gamma function*](https://doi.org/10.5206/mt.v3i1.14591),
Maple Transactions 3(1), Article14591, Eq28 and Theorem3.3/Eq35,
printed pp.13–14. At derivative order1 these give

$$
\psi(z)=\log z-\frac1{2z}
-\sum_{k=1}^{J-1}\frac{B_{2k}}{2kz^{2k}}+R'_J(z),\qquad
|R'_J(z)|\le
\left(1+\sqrt{\pi(J+1/2)}\right)
\frac{|B_{2J}|}{2J|z|^{2J}}
$$

for $J\ge1$, $\Re z\ge0$, $z\ne0$. The primary PDF and these theorem
conditions were inspected. FLINT supplies Bernoulli numbers, finite
power series, Lerch $\Phi$ and Hurwitz $\zeta$ directly. No replacement
algorithm for these functions is implemented.

Set $A=M+1$, $\beta=\pi/\log c$, $a=1/4$, $q=c^{-2}$,
$\delta=c^{-1/2}$, $\alpha=a/\beta$ and
$K=\log c\, (\cosh(\log c/2)-1)/\pi$. For
$P_m=\sum_{r<c}\Lambda(r)r^{-1/2}\sin(2\pi m\log r/\log c)$,
the existing symbol has the exact special-function representation

$$
s_m=-P_m-\frac{Km}{m^2+\alpha^2}
-\frac12\Im\psi(a+i\beta m)
-\frac\delta2\Im\Phi(q,1,a+i\beta m).
$$

The last term is a positive Gamma endpoint correction. The constant
nonoscillatory term is $-\pi/4$. Neither term is dropped.

For $\mu_\ell=\Phi(q,-\ell,a)=\sum_{j\ge0}q^j(j+a)^\ell$,
the finite geometric identity gives

$$
-\frac\delta2\Im\Phi(q,1,a+i\beta m)
=\frac\delta2\sum_{k=0}^{J-1}
\frac{(-1)^k\mu_{2k}}{(\beta m)^{2k+1}}+\rho^E_{J,m},\qquad
|\rho^E_{J,m}|\le\frac{\delta\mu_{2J}}{2(\beta m)^{2J+1}}.
$$

This uses $1/(1+t^2)=\sum_{k<J}(-t^2)^k+(-t^2)^J/(1+t^2)$,
valid for every real $t$; it does not assume the unbounded summation
index $j$ satisfies $(j+a)/(\beta m)<1$. The pole expansion has
remainder at most $K\alpha^{2J}m^{-(2J+1)}$.

To convert the finite digamma sum into inverse-mode coefficients, define

$$
H_J(x)=\log(1+ax/(i\beta))
-\frac{x}{2(i\beta+ax)}
-\sum_{k=1}^{J-1}\frac{B_{2k}}{2k}
\left(\frac{x}{i\beta+ax}\right)^{2k}.
$$

The program requires $\beta>a$, so this function is analytic on a
neighborhood of the closed unit disk. On that disk its modulus is at
most

$$
V_J=-\log(1-a/\beta)+\frac1{2(\beta-a)}
+\sum_{k=1}^{J-1}\frac{|B_{2k}|}{2k(\beta-a)^{2k}}.
$$

Cauchy's coefficient bound controls the Taylor tail by
$V_Jm^{-2J}/(1-1/A)$ for every $m\ge A$. If $[H_J]_{<2J}$
denotes the Taylor polynomial, the program forms

$$
D_J(x)=-\frac\pi4-\frac12\Im_{\rm coeff}[H_J]_{<2J}(x)
-K\sum_{k<J}(-\alpha^2)^kx^{2k+1}
+\frac\delta2\sum_{k<J}(-1)^k\mu_{2k}\beta^{-(2k+1)}x^{2k+1}.
$$

Writing $D_J(1/m)=\sum_j c_j(M/m)^j$, it computes the $c_j$
with FLINT series arithmetic. For $\widetilde s_m=D_J(1/m)-P_m$,
the full common remainder is

$$
|s_m-\widetilde s_m|\le\mathcal A_Jm^{-2J}
+\mathcal B_Jm^{-(2J+1)},
$$

$$
\mathcal A_J=
\frac{(1+\sqrt{\pi(J+1/2)})|B_{2J}|}{4J\beta^{2J}}
+\frac{V_J}{2(1-1/A)},\qquad
\mathcal B_J=K\alpha^{2J}+\frac{\delta\mu_{2J}}{2\beta^{2J+1}}.
$$

The finite order $J$ must pass the output accuracy check; a larger $J$
at fixed $M$ does not guarantee a smaller asymptotic remainder.

For a positive rational $u$, the existing scalar kernel is

$$
\mathsf L_p(u)=M^{p-1}z(u)^A\Phi(z(u),p,A),\qquad
z(u)=e^{2\pi i\log u/\log c}.
$$

For phase1 the program calls Hurwitz $\zeta$ instead. Exact rational
arithmetic merges phases modulo $c^{\mathbb Z}$ and conjugation,
before evaluating any logarithms. In particular $z(3)z(9)=1$ at c27;
its negative sum-frequency terms are retained. If $v_r=\Lambda(r)/\sqrt r$,
the finite centers are

$$
\widetilde{\mathsf T}_p=\sum_jc_j\mathsf L_{p+j}(1)
-\sum_rv_r\Im\mathsf L_p(r),
$$

$$
\begin{aligned}
\widetilde{\mathsf U}_p={}&\sum_{j,k}c_jc_k\mathsf L_{p+j+k}(1)
-2\sum_{j,r}c_jv_r\Im\mathsf L_{p+j}(r)\\
&+\frac12\sum_{r,t}v_rv_t
\left[\Re\mathsf L_p(r/t)-\Re\mathsf L_p(rt)\right].
\end{aligned}
$$

The last sum is ordered. All finite signed products come from the same
symbol; no independent-phase average or unsigned substitute is used.

Set $\epsilon=\mathcal A_J A^{-2J}+\mathcal B_J A^{-(2J+1)}$ and
reuse the existing actual-symbol bound $|s_m|\le B_*$, where
$B_*=2\cosh(\log c/2)+\sum_{r<c}v_r$. Then

$$
F_p=\frac{\mathcal A_J}{M^{2J}}\mathsf L_{p+2J}(1)
+\frac{\mathcal B_J}{M^{2J+1}}\mathsf L_{p+2J+1}(1)
$$

bounds the normalized $\mathsf T_p$ truncation error. The corresponding
$\mathsf U_p$ error is at most $(2B_*+\epsilon)F_p$, since
$|s_m^2-\widetilde s_m^2|\le(2B_*+\epsilon)|s_m-\widetilde s_m|$.
The program adds these balls to all finite-arithmetic and special-function
errors before outward dyadic rounding. All consumers must use the same
moment intervals. Their own coefficient and matrix-norm amplification,
finite middle modes, jet remainder and sign check remain separate.

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/arb_tail_moments.py --power 2 --power 3 --power 4 --power 8 --power 16 --power 32 --power 64 --power 128 --out build/c27-tail-moments.json.zlib
```

Output fields `T`, `U` and `normalized_polynomial` use exact integer
endpoint pairs on the `grid_bits` dyadic grid. The analytic error fields
and `uniform_symbol_error` use `analytic_grid_bits`. Every reported
moment interval has width at most two output grid units, or the program
fails and requests more precision/order. A finite returned special-function
value alone is insufficient. This interface supplies scalar moments;
it supplies no complete coupling Gram, retained-block sign or RH proof.

The saved `tail-moments.json.zlib` packet covers powers
$2,3,4,8,16,32,64,128$, with $M=32768$, $J=32$, 2304 working bits
and a $2^{-640}$ output grid. It uses 3,132 cached scalar kernels;
every one of the sixteen moment intervals has width at most one grid
unit. The included analytic U-errors are all less than $2^{-816}$;
this is only the truncation component of the total exported error.
All six T/U endpoint pairs at powers2,4,128 agree exactly with an
independent 2560-bit execution. The original fourteen-phase sine and
square were checked against the folded expressions at seven modes.
A c3 cutoff shift from64 to71, at powers2,3,8, agrees with the exact
finite original-symbol segment in both moment families. Independent
Hurwitz-progression checks cover production kernels at powers2,4,128.
These checks validate the new interface, without repeating the older
scalar or positivity certificates. Packet SHA-256:
`fc91b487a48f82cb276c0823f28c4a38432f12f822d6e0ed2644a67da7c22ad1`.

## Actual infinite far-coupling Gram

`arb_far_coupling.py` consumes the published actual scalar packet and
the shared normalized moments. It applies the finite geometric expansion
to the actual divided difference, keeping every feature and all signed
products inside each supplied moment. Omitted core entries have an explicit
norm allowance. It returns a rational matrix representation and a uniform
operator error. It does not compute the middle modes or a retained sign.

Work first in the phase-adjusted basis. For $m>M>N$, the original rows are

$$
C_{m0}=-\frac{\sqrt2s_m}{\pi m},\qquad
C_{mn}=\frac{2(ns_n-ms_m)}{\pi(m^2-n^2)},\quad 1\le n\le N.
$$

With $t=M/m$, $x_n=n/M$ and $0\le k<R$, define the exact features

$$
A_k(v)=\frac{v_0}{\sqrt2}\mathbf1_{k=0}
+\sum_{n=1}^N x_n^{2k}v_n,\qquad
B_k(v)=\sum_{n=1}^N s_nx_n^{2k+1}v_n.
$$

The finite geometric jet is

$$
J_m(v)=\frac2{\pi m}\sum_{k=0}^{R-1}t^{2k}
\left(tB_k(v)-s_mA_k(v)\right).
$$

Let $P$ be the feature matrix, ordered as $A_0,\ldots,A_{R-1},
B_0,\ldots,B_{R-1}$. Absolute convergence and finite expansion give

$$
J^*J=\alpha P^*HP,\qquad \alpha=\frac4{\pi^2M},
$$

$$
H_{AA,kl}=\mathsf U_{2k+2l+2},\qquad
H_{AB,kl}=-\mathsf T_{2k+2l+3},\qquad
H_{BB,kl}=\mathsf Z_{2k+2l+4},
$$

where $\mathsf Z_p=M^{p-1}\zeta(p,M+1)$. Thus the full Gram keeps the
same complete symbol in both factors. No cross term is optimized
independently. The jet's zero column is exact at every positive order.
The original centered unphased basis transports this Gram by multiplying
column $n$ of $P$ by $(-1)^n$; the exterior row parity cancels.

The existing finite geometric identity gives, entry by entry,

$$
(C-J)_{m0}=0,\qquad
(C-J)_{mn}=C_{mn}(n/m)^{2R}\quad(n\ge1).
$$

Reuse the existing bound $|s_n|\le B_*$ and put
$q_0=N/M$, $\kappa=1-N/(M+1)$. Since
$|C_{mn}|\le2B_*/(\pi(m-n))$, the integral bound for the positive
inverse-power tail gives

$$
\|C\|^2\le\|C\|_F^2\le
c_2:=\frac{B_*^2}{\pi^2M}\left(\frac{4N}{\kappa^2}+2\right),
$$

$$
\|C-J\|^2\le\|C-J\|_F^2\le
r_2:=\frac{4B_*^2N}{\pi^2\kappa^2M(4R+1)}q_0^{4R}.
$$

Consequently the whole matrix error, rather than just its squared
remainder, is

$$
\|C^*C-J^*J\|\le\eta_{\rm geom}:=2\sqrt{c_2r_2}+r_2.
$$

Only the first $d$ block anti-diagonals need numerical moments. Let
$H^{(d)}$ equal $H$ where $k+l<d$ in each block, and zero elsewhere.
For $q_0>0$, use the exact diagonal scaling
$D_q=\operatorname{diag}((q_0^{2k})_k,(q_0^{2k})_k)$ and put
$\widetilde P=D_q^{-1}P$. The feature definition gives

$$
\|\widetilde P\|^2\le L_R^2:=\frac12+RN(1+B_*^2q_0^2).
$$

The elementary moment bounds
$|\mathsf T_p|\le B_* /(p-1)$,
$0\le\mathsf U_p\le B_*^2/(p-1)$ and
$0<\mathsf Z_p\le1/(p-1)$ follow from the same positive
inverse-power integral bound. Since $B_*\ge1$, every absolute row sum of
$D_q(H-H^{(d)})D_q$ is at most

$$
\tau_d:=\frac{q_0^{2d}}{1-q_0^2}
\left(\frac{B_*^2}{2d+1}+\frac{B_*}{2d+2}\right).
$$

Symmetry gives the same operator-norm bound. Therefore

$$
\|J^*J-\alpha P^*H^{(d)}P\|
\le\eta_{\rm core}:=\alpha L_R^2\tau_d.
$$

When $d=2R-1$, there are no omitted entries and $\tau_d=0$ is used.
For $N=0$ the omitted features vanish and the same zero error applies.
The sparse core need not be positive semidefinite: its missing signed
terms are paid by this scalar allowance. Its first feature row and column
are exact when $d=R$. The default uses $d=R$;
`--full-moments` uses $d=2R-1$.

The input also pays all finite-precision matrix amplification. The exact
rational $P_0$ uses the published scalar midpoints and a dyadic midpoint
for $1/\sqrt2$; its powers $x_n^k$ are exact rationals. Let
$\epsilon_s$ be the maximum original scalar midpoint error and
$\epsilon_a$ the midpoint error for $1/\sqrt2$. Then

$$
\|P\|_F^2\le p_2:=\frac12+
\frac{N(1+B_*^2q_0^2)}{1-q_0^4},\qquad
\|P-P_0\|_F\le e_P:=
\sqrt{\epsilon_a^2+\frac{\epsilon_s^2Nq_0^2}{1-q_0^4}}.
$$

The moment intervals have width at most two $2^{-b}$ grid units. Their
shared exact sparse midpoint matrix $H_0$ therefore satisfies
$\|H^{(d)}-H_0\|\le e_H:=2R\,2^{-b}$. The genuine full Gram $H$ is
positive semidefinite. Masking decreases its Frobenius norm, so its trace
and the original symbol bound give

$$
\|H^{(d)}\|\le\|H^{(d)}\|_F\le\|H\|_F\le H_*:=\sum_{k=0}^{R-1}
\left(\frac{B_*^2}{4k+1}+\frac1{4k+3}\right).
$$

With $p_0=\sqrt{p_2}+e_P$, and a dyadic midpoint $\alpha_0$ with
$|\alpha-\alpha_0|\le\epsilon_\alpha$, the complete rounding allowance is

$$
\eta_{\rm round}=
\alpha\left[H_*(2\sqrt{p_2}e_P+e_P^2)+p_0^2e_H\right]
+\epsilon_\alpha p_0^2(H_*+e_H).
$$

The program evaluates this bound with directed balls. It does not assume
that the midpoint $H_0$ is positive semidefinite. The asserted upper
matrix bound follows after adding the scalar allowance.

For the existing c27 exterior floor, set

$$
h_m=\log(m/\log27)-\frac{\log27}{\pi m}-\frac{59}{10}.
$$

This is the previously supplied full even Neumann–Gamma/prime comparison,
not a new prime estimate. The program requires $h_{M+1}>0$ and rounds
$h_{M+1}^{-1}$ upward to an exact dyadic $D$. Monotonicity yields the
sufficient actual far-debit majorant

$$
\sum_{m>M}h_m^{-1}C_m^*C_m
\preceq D\alpha_0 P_0^*H_0P_0
+D(\eta_{\rm geom}+\eta_{\rm core}+\eta_{\rm round})I.
$$

The scalar allowance controls approximation of $DC^*C$. Replacing the
varying inverse weight by $D$ introduces additional one-sided slack.
Thus this is a contribution to the sufficient diagonal-floor Schur debit;
it does not identify a far block of the full exterior inverse.

The default parameters are $N=1280$, $M=32768$, $R=56$, $b=640$,
with 2304 working bits. The reduced core requires exactly 56 U-moments at
even powers2 through112, 56 T-moments at odd powers3 through113, and 56
Hurwitz moments at even powers4 through114. The seven applicable scalar
moments in the earlier packet are reused. `linear_moment` and
`square_moment` call the same existing supplier formulas separately so
the unused opposite-parity family is not evaluated.

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/arb_far_coupling.py --moment-cache build/c27-far-moments.partial.zlib --out build/c27-far-coupling.json.zlib
```

The program fails if its total weighted operator allowance is not below
$2^{-500}$ or any required moment fails its dyadic output-width check.
`upper_bounds` contains exact nonnegative numerators on the
`bound_grid_bits` grid. `D_numerator` uses `D_grid_bits`; moment and
prefactor interval pairs use `grid_bits`. All pairs use one exact
midpoint with denominator $2^{\text{grid\_bits}+1}$.

The highest-power pair is evaluated first. An optional moment cache saves
completed pairs atomically; reuse requires its explicit `--resume-sha256`
pin and identical source, supplier, seed and numerical parameters. Cached
interval widths and power indices are checked before reuse. The cache
is an intermediate scalar input, not a completed matrix packet. The optional
full core needs 111 moments per family through powers222/223/224.
Its T223 export fails at 2304 and 2560 working bits and passes at 2816;
use `--full-moments --precision 2816` for that case. The default reduced
core does not require this higher-order calculation.

`reconstruct_factors(packet, scalar_packet)` returns exact rational
$\alpha_0,D,P_0,H_0$ and the exported weighted scalar allowance, using
FLINT rational matrices. Its default convention is the centered unphased
basis. This is the explicit far matrix input; a dense 1281-square matrix
and a replacement factorization algorithm are not built here.

The saved `far-coupling.json.zlib` packet contains all 168 required
intervals at the default parameters. Each has width at most one
$2^{-640}$ grid unit. Its exact factors have shapes $112\times1281$
and $112\times112$; seven seed moments are reused unchanged. The
exported weighted allowance has base-two logarithm
approximately $-524.573894274258476$, below the requested $2^{-500}$.
Both highest U112/T113 endpoint pairs agree exactly with independent
2560-bit executions. Validation also covers 768 exact rational row-bound
cases, eight actual finite-vector Gram checks, 32 coordinated input
perturbations, zero-mode, full/reduced selection and pinned-cache guards.
Packet SHA-256:
`ac35edacd7f55de470586e3b8a5915bcaa49b9d2c4c9c2ee073eeaaa8c9c73e7`.

This application reuses the finite geometric mechanism of the existing
Cauchy jets and the actual divided difference. Its new interface is the
actual discrete moment majorant with all matrix errors paid.
The original Fourier/form identification and c27 exterior comparison
remain paper-level inputs. The finite middle rows1281 through32768,
retained sign and cofinal positivity are still missing; no complete
coupling, Lean or RH certificate follows.

## A four-row obstruction to the N1280 diagonal-floor debit

`diagonal_floor_witness.py` reuses the actual scalar table and
`Arithmetic.even_entry`. It supplies only the missing rows1282,1283,1284.
In the original centered unphased cosine basis, take

$$
v=-e_{1277}+2e_{1278}-e_{1280},\qquad \|v\|^2=6.
$$

Write $A_N$ for the actual retained block, $C_m$ for its actual exterior
coupling row, and

$$
d_m=\Re\psi\left(\frac14+\frac{i\pi m}{\log27}\right)
-\log\pi-\frac{59}{10}.
$$

This is the existing full even Neumann–Gamma/prime diagonal lower
comparison; it is sharper than the elementary logarithmic $h_m$ used
in the far input. For each $m=1281,\ldots,1284$, directed evaluation gives

$$
\frac1{16}<d_m<\frac1{15},\qquad 0<h_m<\frac1{15}.
$$

For either choice $w_m=d_m$ or $w_m=h_m$, the sufficient diagonal-floor
Schur target $S_w=A_N-\sum_{m>N}w_m^{-1}C_m^*C_m$ therefore satisfies

$$
\langle v,S_wv\rangle
\le\langle v,A_Nv\rangle
-15\sum_{m=1281}^{1284}|C_mv|^2
<-\frac{13}{5}.
$$

The omitted debit rows are positive semidefinite, so this exclusion
needs no further middle or far summation. The actual retained quadratic
is positive:

$$
\langle v,A_Nv\rangle=30.70294244311003905\ldots,
\quad15\sum_{m=1281}^{1284}|C_mv|^2=33.30890128199787908\ldots.
$$

The saved `diagonal-floor-witness.json` encloses their difference by
$-2.60595883888784003\ldots$ on a $2^{-512}$ grid, using 768 working
bits and the original 640-bit scalar grid. It retains all three newly
generated scalar rows, the four row responses, both floor intervals,
the integer vector and input/source hashes. Result SHA-256:
`f45f3a333b50aa60c38c5b2ee9be416a165a1c719ba7dfb3ca6dfbe251bf1840`.

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/diagonal_floor_witness.py --out build/c27-diagonal-floor-witness.json
```

This witnesses failure of this specified sufficient comparison at N1280,
including the uniform inverse-floor charge16. It is not a negative Weil
test or an RH counterexample. The existing high-complement floor and
far-input majorant remain valid. Further summation alone cannot make
$S_w$ positive; another retained cutoff, a stronger common exterior
comparison or a correlated inverse estimate requires its own analysis.

## Actual finite-column inverse-residual consumer

The existing inverse-residual identity in Liu, Appendix B.4, is applied
to the same actual exterior operator $H=H_{\rm ext}$ after mode1280.
This application assumes the existing actual Fourier/form identification,
the selected Fourier vectors in $D(H)$, and the existing comparison
$H^{-1}\preceq W$. These are analytic inputs, not conclusions of the
numerical programs. The retained space still has1281 coordinates.

Here $W$ is the positive diagonal operator with exact upward dyadic
weights $w_m\ge h_m^{-1}$ for $1281\le m\le32768$, on a256-bit grid.
For $m>32768$ its weight is the already archived constant
$D\ge h_{32769}^{-1}$. The elementary logarithmic $h_m$ is increasing
and positive on this exterior. This reuses the existing comparison;
it does not identify $W$ with the actual exterior inverse.

Let $E$ select the actual exterior columns1281 through1284, put
$Z=HE$, $F=E^*HE$ and $B=E^*Cv$, and reuse the fixed integer vector
$v=-e_{1277}+2e_{1278}-e_{1280}$. The complete joint Gram is

$$
\begin{pmatrix}G&J^*\\J&L\end{pmatrix}
=[Cv,Z]^*W[Cv,Z].
$$

For each prefix of one, two or four selected columns, set $K=L-F$
and $T=J-B$. FLINT proposes a solution of $Kx=T$, which is rounded
to an exact64-bit dyadic vector. The program evaluates the existing
residual bound on this fixed vector:

$$
U(x)=G-2x^*T+x^*Kx
=2x^*B-x^*Fx+(Cv-Zx)^*W(Cv-Zx).
$$

Under the stated analytic inputs,
$\langle Cv,H^{-1}Cv\rangle\le U(x)$. This does not require the
dyadic trial to be an exact optimizer. All signed cross terms remain
in the same actual residual; the selected exterior diagonal uses its
actual scalar formula.

`middle-scalars.json.zlib` adds only the31483 previously missing
scalar pairs1285 through32767, using the unchanged `Arithmetic.scalars`.
The original1281 and32768 rows and the earlier witness rows1282–1284
are reused. Working precision is768 bits and the grid remains640 bits;
all62966 new scalar intervals have width at most one grid unit.
The first100 new rows were reused from a pinned cache. A new cache is
optional; a fresh run computes only missing rows. Packet SHA-256:
`55212d0181673c015de503928265509e594eb0f81b85e3afb34cc28a103d18db`.

Every middle row1281 through32768 contributes to the joint Gram.
Beyond32768 the consumer reuses all168 published signed moment
intervals and the existing112-feature reconstruction. Extending the
feature columns through1284 requires a new call to `far_bounds`;
the former allowance is not transferred. The new weighted far matrix
allowance has $\log_2\epsilon\simeq-524.0636624$.
The physical column combination has squared norms $(6,1,1,1,1)$;
entrywise balls with radii $\epsilon\|l_i\|\|l_j\|$ pay the complete
two-sided far error before applying the dyadic trial.

The saved directed intervals give:

| Exterior trial columns | Upper-debit expression $U(x)$, approximate | Retained minus debit, approximate |
|---|---:|---:|
| none | 50.9913204560 | -20.2883780129 |
| 1281 | 36.3114689536 | -5.6085265105 |
| 1281–1282 | 24.6961199103 | 6.0068225328 |
| 1281–1284 | 18.9893075926 | 11.7136348505 |

These displayed decimals summarize directed intervals; the exact
integer endpoints and dyadic trials are in
`inverse-residual-pilot.json.zlib`, SHA-256
`6db5093bf78fe1dcf43cbbe1b39ef2b1faf4174cf7cf33fe796c0d10ad4d231b`.
Eight distributed new scalar rows agree exactly at1024 bits. An
independent1024-bit five-column calculation, including direct full
residual evaluation, agrees on every scientific output interval,
all31488 weights and the exact trial integers. The serialized analytic
conditions are explicit. Input hashes, source identity, incomplete
coverage, occupied outputs and inadequate precision have rejection guards.

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/arb_middle_scalars.py --out build/c27-middle-scalars.json.zlib
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/arb_inverse_residual_pilot.py --out build/c27-inverse-residual-pilot.json.zlib
```

This is an efficacy test of the estimate on a known obstruction to the
weaker diagonal scheme. Positivity of the actual Schur quadratic on
this high-frequency vector already follows from the existing G7
comparison after1203. It is not a new positivity theorem for that
direction, a sign certificate for the complete retained space, a Lean
verification or an RH proof.

## Full retained-space floating proposals

`float_inverse_residual.py` reuses the same scalar packets, middle
weights and far moments. It extends the actual exterior trial through
1408 and forms the midpoint joint Gram for all1281 retained columns.
NumPy supplies the finite solves and symmetric eigensolver; no new
general matrix algorithm is implemented. For64 and128 exterior columns
it evaluates the optimized midpoint expression
$A-G+T^*K^{-1}T$.

| Exterior columns | Smallest computed eigenvalue | Computed eigenvalues with absolute value at most $10^{-7}$ |
|---|---:|---:|
| 64 | $-3.9827725420\times10^{-6}$ | 23 |
| 128 | $-2.4967957659\times10^{-11}$ | 24 |

The data are in `inverse-residual-float.json`. Neither calculation
pays scalar, far-feature or floating arithmetic errors. The reported
negative values are candidates, not sign certificates; the128-column
case does not establish nonnegativity. The file also retains finite
dyadic candidates for subsequent directed tests. Restricting the
64-column minimum eigenvector to modes0–127 gives a positive computed
quadratic, about $7.64\times10^{-5}$, so that truncated candidate does
not witness its full-space negative value.

```sh
uv run --no-project --python 3.13.12 --with numpy==2.2.6 --with python-flint==0.9.0 --with threadpoolctl==3.6.0 python docs/reports/weil-c27/float_inverse_residual.py --out build/c27-inverse-residual-float.json
```

The directed sections below pay the errors on complete retained
directions and exclude both specified sufficient families. Actual
Schur positivity still requires a different comparison; the successful
fixed-vector pilot does not settle it.

## A directed obstruction to the fixed64-column certificate family

`arb_inverse_residual_obstruction.py` reuses the complete scalar packets,
the archived diagonal weights $W$, all168 far moments and the full1281
coordinates of the64-column candidate. The exact direction and64-entry
trial are in `inverse-residual-k64-candidate.json`. The direction is not
truncated and no full retained Gram or new eigenvector is produced.

Keep $E$ equal to the actual exterior columns1281 through1344. With
$Z=HE$, $F=E^*HE$, $B=E^*C$, put

$$
G=C^*WC,\quad J=Z^*WC,\quad L=Z^*WZ,\quad K=L-F,\quad T=J-B.
$$

For the fixed dyadic direction $v$ and trial $x$, the existing residual
identity gives

$$
m(x)=Q(v-Ex)-\|Cv-Zx\|_W^2
=v^*(A-G)v+2\operatorname{Re}(x^*Tv)-x^*Kx.
$$

The first64 weighted exterior rows alone certify $K\succeq40I$.
Their contribution is $F^*W_{64}F-F$; all omitted rows of $Z^*WZ$
are positive semidefinite. FLINT encloses all64 eigenvalues of this
finite lower matrix with real parts above40. Completing the square
and using the same residual vector $r=Tv-Kx$ yields

$$
\begin{aligned}
v^*(A-G+T^*K^{-1}T)v
&=m(x)+r^*K^{-1}r\\
&\le m(x)+\|r\|^2/40.
\end{aligned}
$$

The directed calculation includes every middle row1281 through32768.
It evaluates the far residual and its64 gradient components using the
existing112-feature core, extended only through1344. The recomputed
operator error $\varepsilon$ pays $\varepsilon\|v-Ex\|^2$ for the
quadratic and $\varepsilon\|v-Ex\|$ for each gradient coordinate.
The gradient is $Z^*W(Cv-Zx)-Bv+Fx=Tv-Kx$; the selected residual
rows are included in this subtraction.

The768-bit execution, with outward512-bit output intervals, gives

| Quantity | Approximate enclosed value |
|---|---:|
| Fixed trial margin $m(x)$ | $-3.98277254117260544336127915956\times10^{-6}$ |
| Residual correction $\|Tv-Kx\|^2/40$ | $1.671986500522932\times10^{-32}$ |
| Upper bound for the optimized-family quadratic | $-3.98277254117260544336127914284\times10^{-6}$ |

The last interval has strictly negative upper endpoint. Thus **every**
matrix trial $X$ for this fixed $E,W$ fails to make $A-U(X)$ positive
semidefinite: its quadratic on this same $v$ is at most the displayed
optimized value. Improving the finite solve or changing coordinates
cannot repair this specified family. A different $E$ or inverse
majorant $W$ remains a separate task.

The original complete Fourier/form identification, exterior domain and
$H^{-1}\preceq W$ remain analytic premises. This excludes a sufficient
lower-certificate family; it neither gives negative actual Weil energy
nor determines the sign for128 exterior columns or RH. The experiment
uses existing square completion, inverse-residual comparison and FLINT
algorithms, without a new general theorem, originality or Lean claim.

The packet is `inverse-residual-k64-obstruction.json`, SHA-256
`d2a447699a7893be537c461bab3d43b619d0611d7dc90052fb9b9df571d4ca73`.

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/arb_inverse_residual_obstruction.py --out build/c27-k64-obstruction.json
```

The optional `--propose` mode supplies only the64 trial coefficients for
the fixed direction using NumPy. It requires the same pinned NumPy and
threadpoolctl versions as the floating supplier above. Verification
uses the saved exact trial and requires neither NumPy nor another
full-space solve.

## The128-column family has a directed obstruction too

The same consumer accepts the separate128-column candidate, with actual
exterior modes1281 through1408. Only the selected-column count changes;
the complete1281-coordinate direction is reused from that existing
proposal. The64-column packet and its executed producer remain pinned
to their original commit. No old scalar, moment, Gram or eigensolver
output is regenerated.

For this instance the first128 weighted rows again certify $K\succeq40I$.
The same square-completion argument above gives:

| Quantity | Approximate enclosed value |
|---|---:|
| Fixed trial margin | $-2.49668870337129250789512679741\times10^{-11}$ |
| Residual correction | $7.169494010542168\times10^{-34}$ |
| Optimized-family quadratic upper | $-2.49668870337129250789505510247\times10^{-11}$ |

The upper endpoint is strictly below $-2.49\times10^{-11}$, excluding
every $X$ for this specified128-column $E$ and the same $W$. This is a
separate complete direction and trial; the64-column sign was not used
to infer this sign. Actual Schur/Weil positivity and RH remain unresolved.

Extending the existing far features through1408 requires operator error
$\varepsilon\simeq6.0325461\times10^{-154}$. Its quadratic interval
has width about16.17664 units on the $2^{-512}$ grid, exceeding that
exporter's two-unit contract. On the $2^{-504}$ grid its width is
about0.06319 units. The consumer keeps its full error allowance and
strict raw-ball negative test; this instance uses504-bit exported
intervals with768 working bits.

Exact inputs and result are `inverse-residual-k128-candidate.json` and
`inverse-residual-k128-obstruction.json`. The result SHA-256 is
`3b850c5e42f0c5671884a2abdabbc66b5e7611df9bfb6e24f63b512215679a51`.

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/arb_inverse_residual_obstruction.py --candidate docs/reports/weil-c27/inverse-residual-k128-candidate.json --output-bits 504 --out build/c27-k128-obstruction.json
```
