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

## An actual inverse-residual correction beyond the excluded span

Keep the exact full dyadic direction $v$ and trial $y_0=E_{128}x_0$
in `inverse-residual-k128-candidate.json`, the same actual exterior $H$
and the archived $W\succeq H^{-1}$. Let $V$ select exterior modes1281
through2048. The existing dense linear solve proposes a correction from

$$
(V^*HV)z=V^*(Cv-Hy_0),\qquad y=y_0+Vz.
$$

This is an actual-inverse proposal, rather than another optimization of
the excluded $Kx=Tv$ system. The saved $y$ has768 exact96-bit dyadic
coordinates, including640 nonzero coordinates beyond the old span.
The finite solve supplies no infinite-inverse certificate by itself.

`inverse_residual_enrichment.py` accepts the saved finite trial through
the existing full-exterior residual inequality. Put $r=Cv-Hy$ and
$\ell=Q(v-y)-\|r\|_W^2$. Under the same common-form/domain and inverse
comparison premises,

$$
\ell\le v^*(A-C^*H^{-1}C)v\le Q(v-y).
$$

The consumer reuses the old directed $Q(v-y_0)$ and evaluates its exact
correction, with $\delta=y-y_0$:

$$
Q(v-y)=Q(v-y_0)-2\operatorname{Re}\langle\delta,Cv-Hy_0\rangle
+\langle\delta,H\delta\rangle.
$$

Every middle row1281–32768 is included once; the first768 rows serve
both this correction and the new residual. All168 old far moments are
reused. The112-feature columns are extended through2048 with the
complete symbol bound $2H_{\rm gamma}+\sum w_{p^j}$ and freshly paid
operator allowance approximately $1.80743\times10^{-135}$. No old
scalar, moment, retained Gram or eigensolver production is repeated.

| Quantity on the same full dyadic $v$ | Approximate enclosed value |
|---|---:|
| Old complete residual energy | $4.3635328423409634\times10^{-10}$ |
| New complete residual energy | $1.0904887215065086\times10^{-10}$ |
| Old sufficient lower target | $-2.4966887033712925\times10^{-11}$ |
| New sufficient lower target $\ell$ | $2.8606682320961359\times10^{-10}$ |
| New upper trial energy $Q(v-y)$ | $3.9511569536026445\times10^{-10}$ |

Exact endpoints give new residual energy below one quarter of the old
energy and the conditional scalar enclosure
$2.86\times10^{-10}<v^*(A-C^*H^{-1}C)v<3.96\times10^{-10}$.
This repairs the former witness direction after changing the exterior
trial. It does not undo either fixed-family exclusion or certify the
other retained directions, the full1281-dimensional matrix, cofinal
support positivity or RH. The vector has its original exact mass;
these quadratics are not claimed to be normalized eigenvalues.
No new inverse algorithm, general theorem or Lean verification is claimed.

The candidate and result are
`inverse-residual-enrichment-2048-candidate.json` and
`inverse-residual-enrichment-2048-result.json`. Result SHA-256:
`1084eca5cb36a75a5a6d82fe9cc4fa2ef88ab0b59b179d7b7fa29729f7aace22`.
The result uses768 working bits and384-bit exported intervals;
the far allowance is separately exported on a640-bit grid.

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/inverse_residual_enrichment.py --trial docs/reports/weil-c27/inverse-residual-enrichment-2048-candidate.json --out build/c27-inverse-enrichment.json
```

Without `--trial`, the program proposes the finite correction using
NumPy2.2.6 and threadpoolctl3.6.0 in addition to the pinned python-flint.
The directed consumer needs neither NumPy nor another solve.

## A complete joint inverse trial and its reduced obstruction

The finite solver interface reuses Poloni's Trummer-like displacement
machinery, [arXiv:0903.4569v3](https://arxiv.org/pdf/0903.4569v3):
the same-node definition and finite Schur update on p.3, Theorem4 on
p.11, the separate diagonal update on p.12 and Algorithms5–6 on
pp.21–22. The actual diagonal is additional input; generators alone do
not recover it. Section7, p.17, warns about generator growth and actual
residual checks. No source implementation is copied, general algorithm
invented or universal stability inferred.

`structured_inverse_residual_interval.py` uses all1281 original
orthonormal retained coordinates, the same actual exterior $H$ and
$W\succeq H^{-1}$. Select exterior modes1281–2048 by $V$, and set
$F=V^*HV$, $B=V^*C$, $X=F^{-1}B$, $Y=VX$. The trial is this exact
implicit finite inverse, rather than the previous96-bit rounded scalar
trial. Two right-hand sides $F^{-1}G_E$ reconstruct every column of $X$
from the finite rank2 displacement identity. The selected768 residual
rows vanish exactly. The remaining finite rows2049–32768 are all kept;
their weighted residual Gram has displacement rank at most4, with
actual diagonals accumulated separately. The finite product identity
has the plus sign

$$
D_P(R^*WR)-(R^*WR)D_P
=\frac2\pi\bigl(TJ_2(G'_P)^*+G'_PJ_2T^*\bigr),\qquad
T=R^*WG'_O.
$$

All168 archived far moments and112 features remain. The transformed
features are $P_P-P_EX$. The full far allowance, approximately
$1.80743\times10^{-135}$ before transport, pays
$\varepsilon(1+\|X\|_F^2)<2.07733\times10^{-135}$.
The complete sufficient target uses228 finite generator columns and
the separately supplied diagonal. No infinite generator sum, bounded
infinite displacement or positive signed-moment midpoint is assumed.

The reusable640-grid target is `structured-joint-2048-target.json.zlib`,
SHA-256 `ea4a761f76b59a11b2d04faf80531901a5164f4145c2424a4e79889ae3f72e22`.
Alternative factorization consumers can use this target directly.

| Saved result | Conclusion |
|---|---|
| `structured-joint-2048-float.json` | Smallest reading about $-4.83\times10^{-12}$; finite-middle subset discrepancy about $4.62\times10^{-12}$, so no sign or small-eigenspace claim. |
| `structured-joint-2048-pivots.json` |1263 positive pivots, then uncertain mode24;18 coordinates remain. |
| `structured-joint-2048-dense-schur.json` | Validated dense solve of the positive block; eigenvalue isolation inconclusive, but all18 reduced diagonal upper endpoints negative. |
| `structured-joint-2048-diagonal-check.json` | Exact integer consumption of the same stored endpoints rejects this explicit sufficient lower target. |

The first remaining coordinate, mode0, has reduced diagonal upper
endpoint less than $-5.77\times10^{-31}$. A negative Schur diagonal
already precludes positive semidefiniteness of this target because the
eliminated block is positive. This is not a claim of18 negative
eigenvalues, negative actual $A-C^*H^{-1}C$, or failure of every trial
on the selected support. In particular the target includes its paid
conservative far allowance; rejection does not remove that allowance.
The common-form, domain and actual infinite-inverse comparison remain
the same analytic premises. Full c27 positivity, cofinal positivity
and RH are unresolved; no Lean verification is claimed.

The rational `--self-check` covers the finite residual/Gram identities,
signed low-rank update and dense LDL correspondence. A separate
1024-bit run reproduced every target endpoint and the same1263-pivot
prefix; its324 reduced intervals overlap the saved768-bit enclosures.
This is precision reproduction by the same serial native reviewer,
not independent algorithm or model diversity. No old scalar, moment,
Gram, obstruction or enrichment production is repeated.

The stored target and reduced enclosures suffice for the endpoint
conclusion, using only the standard library:

```sh
python3 docs/reports/weil-c27/structured_schur_diagonal_check.py --target docs/reports/weil-c27/structured-joint-2048-target.json.zlib --pivots docs/reports/weil-c27/structured-joint-2048-pivots.json --reduced docs/reports/weil-c27/structured-joint-2048-dense-schur.json --out build/c27-joint-diagonal-check.json
```

For a new production or a different factorization, use fresh outputs:

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/structured_inverse_residual_interval.py --self-check
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/structured_inverse_residual_interval.py --target-out build/c27-joint-target.json.zlib --out build/c27-joint-pivots.json
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/structured_inverse_residual_interval.py --target docs/reports/weil-c27/structured-joint-2048-target.json.zlib --out build/c27-joint-reused-pivots.json
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/structured_residual_schur_finish.py --target docs/reports/weil-c27/structured-joint-2048-target.json.zlib --pivots docs/reports/weil-c27/structured-joint-2048-pivots.json --out build/c27-joint-dense-schur.json
```

## Optimizing the768-column support still leaves an obstruction

The preceding negative diagonal excludes one chosen trial, not all
trials on that support. `structured_support_test.py` addresses the
missing optimization using the existing square-completion argument
already used for the64- and128-column families.

Let $L_0$ be the saved complete target, with its1263-coordinate positive
principal block $P$. Lift unresolved Schur coordinate0 as
$(-P^{-1}be_0,e_0)$, then round the retained coefficients to exact
256-bit dyadics. The complete1281-coordinate vector $v$ is stored in
`support-2048-direction.json`; its exact mass is retained, with
$\|v\|^2\simeq2.0458788430940653$. Its saved-target quadratic is strictly
below $-5.77\times10^{-31}$. It is a new direction, distinct from the
earlier128-column witness which the scalar enrichment repaired.

Keep $V$ selecting modes1281–2048 and the same actual $H,C,W$. Put

$$
Z=HV,\quad F=V^*HV,\quad
K=Z^*WZ-F,\quad t=Z^*WCv-V^*Cv.
$$

For $m(x)=Q(v-Vx)-\|Cv-Zx\|_W^2$ the existing identity gives

$$
\sup_xm(x)=m(x_0)+(t-Kx_0)^*K^{-1}(t-Kx_0)
\le m(x_0)+16\|t-Kx_0\|^2.
$$

Here $K\succeq I/16$ is certified by768 positive pivots of the finite
lower matrix formed using **all** middle rows1281–32768; the omitted
actual far Gram is positive semidefinite. No positivity of the signed
moment midpoint is used. Finite Gram generators use the same product
rule as above, with selected rows overlapping the column set and
their actual diagonals supplied separately.

The expensive finite sums are archived in
`support-2048-finite-action.json.zlib`. Alternative floor/optimization
consumers use `--finite` and rescan no middle rows. The far consumer
reuses all168 moments, the112-feature core and the already paid
extended-support error from the original512-grid pivot packet.
It pays $\varepsilon\|v-Vx_0\|^2$ for the quadratic and
$\varepsilon\|v-Vx_0\|$ for each gradient coordinate. A finite solve
only proposes $x_0$; the stored trial has768 exact384-bit dyadic
coefficients, and the full residual comparison accepts the result.

| Quantity on this same $v$ | Approximate enclosed value |
|---|---:|
| Fixed trial margin $m(x_0)$ | $-4.22459205900814696498591\times10^{-31}$ |
| Optimizer correction bound on the stored grid | $\le2^{-512}<7.46\times10^{-155}$ |
| Upper bound for $\sup_xm(x)$ | $-4.22459205900814696498591\times10^{-31}$ |

The exact upper endpoint is below $-4.22\times10^{-31}$. Thus every
trial matrix on this fixed768-column support, with this same $W$,
fails to make its inverse-residual lower target positive semidefinite.
This strengthens the preceding chosen-trial rejection; it neither
gives negative actual Weil energy nor excludes other supports or
inverse majorants. The same common-form/domain/inverse-comparison
premises remain. Full positivity, cofinal positivity and RH are open;
these quadratics are not normalized eigenvalues or Lean results.

The result is `support-2048-optimized.json`. The new exact rational
check includes overlapping row/column sets and the signed $K$-floor
LDL correspondence. The saved direction, finite action and result
allow later consumers to reuse all numerical production.

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/structured_support_test.py --self-check
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/structured_support_test.py --candidate docs/reports/weil-c27/support-2048-direction.json --finite docs/reports/weil-c27/support-2048-finite-action.json.zlib --out build/c27-support-reused.json
```

Fresh direction/finite production, if required, uses `--prepare` and
`--finite-out`; it does not regenerate old scientific packets:

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/structured_support_test.py --prepare --out build/c27-support-direction.json
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/structured_support_test.py --candidate build/c27-support-direction.json --finite-out build/c27-support-finite.json.zlib --out build/c27-support-optimized.json
```

## Retain positive even Gamma boundary channels

`gamma_boundary_update.py` changes the inverse majorant while keeping
exactly the saved complete direction and saved768-coordinate trial.
It reuses the existing even Gamma boundary decomposition, variational
inversion, Woodbury identity and inverse-residual comparison. No new
general inverse algorithm is supplied.

For $d=\log27$, $b_j=2j+1/2$ and $n>1280$, the original orthonormal
cosine coefficients give the positive boundary term

$$
B_j=(1-e^{-b_jd})h_jh_j^*,\qquad
h_{j,n}=\frac{2b_j(-1)^n}{\sqrt d\,[b_j^2+(2\pi n/d)^2]}.
$$

For a finite bin $[a,b)$, set $c_a=1-e^{-b_ad}$ and
$U_{[a,b)}=\sqrt{c_a/(b-a)}\sum_{j=a}^{b-1}h_j$.
Finite Cauchy--Schwarz gives $\sum_{j=a}^{b-1}B_j\succeq
U_{[a,b)}U_{[a,b)}^*$. Nineteen disjoint bins use the endpoints
$0,1,2,3,5,8,13,21,34,55,89,144,233,377,610,987,1597,2584,4181,6765$.
The Fibonacci endpoints are a choice of scale grouping; their necessity
and optimality are not asserted.

The archived upward reciprocal weights imply
$H\succeq W^{-1}+UU^*$ on the same actual even form domain. Existing
variational inversion and Woodbury therefore yield

$$
H^{-1}\preceq W_{\rm G}:=W-WU(I+U^*WU)^{-1}U^*W\preceq W.
$$

These columns have infinite Fourier support. This changes the fixed
$W$ in the preceding exclusion, while retaining its $v,x_0$ and actual
$H$. For $r=Cv-HVx_0$, $\beta=U^*Wr$ and $S=I+U^*WU$, the improvement
to that same trial margin is $\beta^*S^{-1}\beta$.

`support-2048-residual-action.json.zlib` archives the previously unsaved
actual residual on all31488 middle rows. The original scalar tables,
exact trial, complete margin and far residual energy are consumed.
Rows through4096 are direct; subsequent rows use the existing finite
geometric expansion of $1/(m^2-n^2)$ to order128, including the
absolute-coefficient remainder. The two coefficient polynomials use
the signed combinations $(-1)^n n s_nu_n$ and $(-1)^nk_nu_n$, with
$u=(v,-x_0)$ and $k_0=1/\sqrt2$, $k_n=1$ otherwise. At
$t=(2048/m)^2$ the remainder radius is

$$
\frac{2}{\pi m^2}\frac{t^{128}}{1-t}
\left(\sum_n|n s_nu_n|+m|s_m|\sum_n|k_nu_n|\right).
$$

Four direct rows check the geometric evaluator at4097,8192,16384,32768;
the complete finite energy overlaps its independent archived value.
The residual preparation took8.78 seconds. It reruns no old Gram,
target, solve, scalar or far-moment producer.

For the new thin projections, FLINT digamma evaluates the finite rate
sums. Each column satisfies $|U_{k,m}|\le a_k/m^2$, where
$a_k=2\sqrt{c_a/((b-a)d)}(b-a)(a+b-1/2)(d/(2\pi))^2$.
With archived constant far weight $D$ and actual far energy upper bound
$E_\infty$, put $\eta_k=Da_k^2/(3M^3)$ and
$\delta=\sum_k\eta_k$. Weighted Cauchy gives mixed-tail error
$\sqrt{\eta_kE_\infty}$ per coordinate, and the positive far Gram
is bounded by $\delta I$ using its trace. Entrywise bounds are not
promoted to rank-one Loewner bounds.

Write $S_0=I+U_0^*W_0U_0$, $\beta_0=U_0^*W_0r_0$ for all middle
rows. An exact stored256-bit dyadic $z$ gives

$$
\beta^*S^{-1}\beta\ge2\Re z^*\beta_0-z^*S_0z
-2\sum_k|z_k|\sqrt{\eta_kE_\infty}-\delta\|z\|^2.
$$

An upper enclosure follows from $S\succeq S_0\succeq I$ and the
Euclidean norm of the mixed-tail error. The512-bit directed consumer,
with256-bit output endpoints, took7.68 seconds and gives:

| Quantity on the unchanged complete direction and trial | Enclosed value |
|---|---:|
| Original margin | $-4.224592059008147\times10^{-31}$ |
| Lower gain | $6.603392989188944\times10^{-31}$ |
| Upper gain | $6.703543517093231\times10^{-31}$ |
| Updated margin | between $2.378800930180797\times10^{-31}$ and $2.478951458085084\times10^{-31}$ |

The exact lower endpoint exceeds $2.37\times10^{-31}$, so the actual
retained Schur quadratic is strictly positive on this same direction
under the original common-form/domain identifications, existing even
boundary comparison and sound reused enclosures. This certifies one
direction; it does not optimize the new trial family or certify the
whole matrix. Full c27 positivity, cofinal positivity and RH remain
unresolved. This is a numerical application, without Lean verification
or a literature-priority claim.

The new residual is reusable by later thin consumers:

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/gamma_boundary_update.py --residual docs/reports/weil-c27/support-2048-residual-action.json.zlib --out build/c27-gamma-boundary.json
```

Only when a fresh residual archive is required:

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/gamma_boundary_update.py --prepare --out build/c27-residual-action.json.zlib
```

## A common Gamma update for every retained coefficient

`joint_gamma_boundary.py` applies the existing positive Gamma channels
and inverse-residual comparison to the original common exact trial
$Y=EF^{-1}B$, not the saved scalar optimizer of the preceding section.
The selected exterior modes are1281–2048, while all original retained
modes0–1280 and every actual exterior row are preserved. The scalar
rows, diagonal weights, complete old structured target,112-feature far
covariance and its error bounds are reused.

Write $R=C-HY$, $Q_Y=A-C^*Y-Y^*C+Y^*HY$ and
$T_{\rm G}=Q_Y-R^*W_{\rm G}R\preceq S_N$. This ideal sufficient
target is distinct from the actual Schur matrix. With the row split
at32768, the new common projections are

$$
A_0=U_0^*W_0U_0,\qquad \Beta_0=U_0^*W_0R_0.
$$

`joint-residual-generators-2048.json.zlib` retains the previously
unsaved two-generator residual action on rows2049–32768. It reuses the
archived selected block, solves only two generator right sides, and
uses direct rows through4096 and a384-term finite geometric expansion
thereafter. Its absolute remainder and four direct geometric samples
are included. `joint-gamma-6765-finite.json.zlib` retains the complete
$19\times19$ Gram and $19\times1281$ thin projection. All31488 middle
rows enter $A_0$, including the selected rows. The selected-row
residual is exactly zero, so $\Beta_0$ uses the30720 outside rows.
The512-term moment expansion has paid uniform remainder with ratio
$1280/2049$. No old scalar or far-moment production is repeated.

The actual far residual Gram has the existing complete majorizer
$G_R=P^*KP+aI$. Only this sum, not the signed core alone, is used as
an upper bound for the positive Gram. The Gamma far trace bound gives
$U_\infty^*W_\infty U_\infty\preceq\delta I$, with
$\delta\le0.000106746689$. For a common $\rho>0$, finite Cauchy and
the existing variational inverse identity yield

$$
L_\rho=L_{\rm old}
+\Beta_0^*[I+A_0+(1+\rho^{-1})\delta I]^{-1}\Beta_0
-\rho G_R\preceq T_{\rm G}\preceq S_N.
$$

The factor $1+\rho^{-1}$ pays both the mixed tail and the negative
far-channel Gram. The signed structured kernel becomes
$\operatorname{diag}((1+\rho)K,-[I+A_0+(1+\rho^{-1})\delta I]^{-1})$;
the original four displacement generators and independent diagonal
are preserved. This applies existing identities and algorithms to
new thin data; it introduces no general inverse algorithm or
literature-priority claim.

The original target's matrix endpoints use a640-bit grid. Its
allowance metadata uses the512-bit grid explicitly supplied by the
original pivot packet. The new consumer validates the inherited
metadata, decodes the two grids separately, and judges the same
outward target it archives. The corrected `rho32-v2` packets below
are the retained derived outputs.

At $\rho=1/32$,1263 positive pivots are certified and the next interval
contains zero at mode24. The existing validated dense finish leaves
18 coordinates. Eigenvalue isolation is inconclusive; all18 reduced
diagonal upper endpoints are nevertheless negative. One such
endpoint rejects this lower target, without counting negative
eigenvalues or determining actual Weil positivity.

`structured_residual_schur_finish.py --actions-out PATH` saves its
already computed positive-block inverse action. The consumer
`joint_gamma_direction.py` uses that action at remaining mode0 to
build a full exact256-bit dyadic direction, with squared norm
approximately2.058687316. It repeats no dense solve or row production.
For $\beta_0=\Beta_0v$, $S_0=I+A_0$, $E_v=v^*G_Rv$ and
$e=\sqrt{\delta E_v}$, its independent gain bounds are

$$
2\Re z^*\beta_0-z^*S_0z-2e\|z\|-\delta\|z\|^2
\le\beta^*S^{-1}\beta
\le\beta_0^*S_0^{-1}\beta_0+2e\|\beta_0\|+e^2.
$$

The complete old target differs from the ideal old sufficient target
by between zero and $2a\|v\|^2$. With that allowance included, the
exact directed endpoints establish

$$
-1.18\times10^{-31}<v^*T_{\rm G}v<-1.14\times10^{-31}.
$$

The common-Young target on the same dyadic direction is approximately
$-1.177716569\times10^{-31}$. Thus changing only $\rho$ or sharpening
the tail payment cannot make this fixed $Y,W_{\rm G}$ target positive.
This excludes the ideal sufficient target for this common trial,
not all trials on the same selected support, additional Gamma
channels, a sharper inverse majorant, the actual Schur matrix or RH.
The preceding positive scalar repair uses a different direction and
trial and remains valid. These results are conditional numerical
applications under the existing actual form/domain and comparison
premises, without Lean verification. Full c27 positivity and cofinal
positivity remain unresolved.

The retained inputs and consumers permit reuse without regeneration:

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/joint_gamma_boundary.py --finite docs/reports/weil-c27/joint-gamma-6765-finite.json.zlib --target-out build/c27-joint-gamma-target.json.zlib --out build/c27-joint-gamma-pivots.json
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/joint_gamma_direction.py --lift-target docs/reports/weil-c27/joint-gamma-6765-rho32-v2-target.json.zlib --lift-pivots docs/reports/weil-c27/joint-gamma-6765-rho32-v2-pivots.json --actions docs/reports/weil-c27/joint-gamma-6765-rho32-v2-actions.json.zlib --out build/c27-joint-gamma-direction.json
```

Only when new channels require new thin projections, reuse the saved
generators with `--produce --generator PATH`; do not repeat the
unchanged residual-generator production. The two production inputs,
corrected target/pivot/Schur packets, saved lift and dyadic-direction
result retain their source hashes and parameter provenance.

## Every trial on the same support with the changed Gamma comparison

`gamma_support_optimization.py` tests the same complete saved dyadic
direction against every trial on exterior modes1281–2048 with the
changed nineteen-channel $W_{\rm G}$. It reuses the old finite support
Gram, original complete target, residual generators, scalar rows and
signed infinite far core. It produces only the missing $19\times768$
finite projection $D_0=U_0^*W_0Z_0$ and two mixed inverse generators,
where $Z=HE$. `gamma-support-6765-finite.json.zlib` archives these
direction-independent actions for later consumers.

Every middle row1281–32768 enters $D_0$, including the selected diagonal.
Rows through4096 are direct. The384-term finite geometric expansion
on the remaining rows has uniform ratio $2048/4097$; its absolute
coefficient remainder is included. The cache's producer source hash
is historical provenance; the subsequent consumer-only changes to
the paid allowance and finite floor do not change its producer formulas.

Let $Y_0=EF^{-1}B$, $F=E^*HE$, $R=C-HY_0$, $G_Z=Z_0^*W_0Z_0$,
$V=F^{-1}G_E$, $\widetilde G_O=G_O-Z_0V$,
$T_Z=Z_0^*W_0G_O$ and $T_R=R_0^*W_0\widetilde G_O$.
The old $G_Z,V,T_Z,T_R$ and transformed retained generator $G_P$ are
already saved. Only $V_A=F^{-1}(T_Z-G_ZV)$ needs two new right sides.
The existing finite displacement identity gives

$$
(F^{-1}Z_0^*W_0R_0)_{mn}=\frac{\kappa}{m^2-n^2}
\left[(V_A)_{m0}(G_P)_{n1}-(V_A)_{m1}(G_P)_{n0}
 +V_{m0}(T_R)_{n1}-V_{m1}(T_R)_{n0}\right].
$$

Selected $m$ and retained $n$ are disjoint; no overlapping-node
diagonal is inferred from this formula. The actual independent
diagonals in $F$ and $G_Z$ remain intact.

For $S_0=I+U_0^*W_0U_0$, dropping the complete nonnegative far
square in the existing channel variational formula gives

$$
K_{\rm G}=Z^*W_{\rm G}Z-F
\succeq K_{{\rm G},0}:=G_Z-F-D_0^*S_0^{-1}D_0.
$$

This keeps the joint channel dependence. The existing structured
interval LDL certifies all768 positive pivots of
$K_{{\rm G},0}-I/16$, using the rank19 update and independent
diagonal. The old fixed-$W$ floor alone is not transferred through
the negative channel update. The result packet saves the complete
new pivot order and intervals.

An exact384-bit dyadic correction $t_0$ is proposed by the changed
finite normal system. For the unchanged complete retained $v$, put

$$
m_{\rm G}(t)=v^*Q_{Y_0}v+t^*Ft-\|Rv-Zt\|_{W_{\rm G}}^2,
\qquad g_{\rm G}(t)=Z^*W_{\rm G}Rv-K_{\rm G}t.
$$

The old paid far allowance is added back as its full interval, on
the original512-bit metadata grid; the matrix grid remains640 bits.
The new far error is paid on the transformed original-column
coefficient norm. The preceding Gamma gain bounds enclose the full
ideal margin. For the gradient, the same Cauchy/resolvent estimates
give $\|z_{\rm true}-z_0\|\le e+\delta\|z_0\|$ and
$\|D_\infty\|\le\sqrt{\delta\operatorname{Tr}G_{Z,\infty}}$.
The trace includes the signed far-core error. Both the finite $D_0$
amplification and full far mixed action are paid; the archived
extended-support error also pays $\sqrt{768}\varepsilon\|u\|$.

The exact directed endpoints in `gamma-support-6765-optimum.json`
give the rounded bounds

$$
-4.102\times10^{-32}<m_{\rm G}(t_0)<-3.652\times10^{-32},
\qquad\|g_{\rm G}(t_0)\|<1.943\times10^{-17}.
$$

Existing square completion with the newly certified $1/16$ floor
then yields

$$
\sup_t m_{\rm G}(t)
\le m_{\rm G}(t_0)+16\|g_{\rm G}(t_0)\|^2
<-3.04\times10^{-32}.
$$

Every trial on this same768-mode support therefore fails to certify
this direction for this fixed $W_{\rm G}$, conditional on the existing
actual common-form, selected-domain, comparison and sound enclosure
premises. This includes every common linear trial on the same
support. The optimized upper is a bound on the best sufficient
target in this family, not a negative upper for the actual Schur
form. Larger supports, additional channels and stronger inverse
comparisons remain open. There is no new generic algorithm,
literature-priority claim, Lean verification, full c27/cofinal
positivity or RH proof.

The cached consumer avoids all old scientific production:

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/gamma_support_optimization.py --finite docs/reports/weil-c27/gamma-support-6765-finite.json.zlib --trial docs/reports/weil-c27/gamma-support-6765-optimum.json --out build/c27-gamma-support-check.json
```

`--trial` reuses the saved exact correction and revalidates the floor,
margin and gradient without another proposal solve. Omit it only
when a new proposal is needed. `--produce --out PATH` is reserved
for a genuinely missing thin cache. The retained cache SHA256 is
`63b908fee9a0ed4f3dfc3e5480d57fcb1678c0c03877c708c027fa1d1cbdff6f`.

## Retain the even pole and the within-bin rate contrasts

`--keep-even-pole` uses the existing even pole without a new thin
matrix. For the singleton rate0, the exterior cosine moment is
$\sqrt2\sinh(d/4)h_{0,n}$. The actual pole
$2|\int f\cosh(x/2)|^2$ therefore adds
$4\sinh^2(d/4)h_0h_0^*$, and

$$
c_0+4\sinh^2(d/4)=e^{d/2}c_0,\qquad
c_0=1-e^{-d/2},\quad d=\log27.
$$

Only the existing channel0 is scaled by $e^{d/4}$: its Gram row and
column, Beta row and support-action row. Its far trace contribution
increases by $(e^{d/2}-1)\eta_0$, where
$\eta_0=Dc_0(d/(2\pi))^4/(3M^3d)$.
The actual $H$, residual and target already contain the pole; this
changes only their inverse comparison. The identity applies on even
tests; the unrestricted pole pair is not asserted positive.

`gamma_channel_resolution.py` applies classical finite Bessel/Haar
projection to the rate bins. For a root bin $[a,b)$ split at $m$,
write $p=m-a$, $q=b-m$, $\ell=p+q$ and retain

$$
V=\sqrt{c_a}\left[
\sqrt{\frac q{p\ell}}\sum_{j=a}^{m-1}h_j
-\sqrt{\frac p{q\ell}}\sum_{j=m}^{b-1}h_j\right].
$$

Use the same original $c_a$ for every detail inside this root.
The old mean and retained details are orthogonal coefficient
projections, so their joint positive matrix remains below the actual
Gamma boundary contribution. At depth1 there are16 new contrasts
and35 channels in total. The new cache contains only their
$19\times16$ cross, $16\times16$ Gram, $16\times1281$ Beta and
$16\times768$ support blocks. Parent aggregates and all scientific
scalar/residual/far suppliers are consumed without regeneration.

The finite Beta action uses the saved transformed outside generator
on rows2049–32768 and512 geometric terms with ratio1280/2049.
The support action uses direct rows1281–4096, including the actual
independent diagonal, and384 terms thereafter with ratio2048/4097.
The absolute finite remainders are paid. For each detail the existing
tail estimate follows from

$$
|V_n|\le\frac{a_V}{n^2},\quad
a_V=2\sqrt{\frac{c_a}{d}}\sqrt{pq\ell}
\left(\frac d{2\pi}\right)^2,\qquad
\eta_V=\frac{Da_V^2}{3M^3}.
$$

The derivative of $b/(b^2+\omega^2)$ has magnitude at most
$1/\omega^2$, and the two rate means differ by $\ell$.
The enlarged far Gram is bounded by its trace, adding the old and
new channel energy allowances while preserving the cross terms.
The consumer assembles the complete common blocks, revalidates
the finite shorting floor and keeps refinement and pole flags in
the exact trial-family provenance. No new generic algorithm or
literature-priority claim is made.

The saved results `gamma-support-pole-6765-optimum.json` and
`gamma-support-pole-contrasts-depth1-optimum.json` each certify all768
positive pivots of the changed finite $K_0-I/16$. The directed
square-completion upper bounds are, respectively,

$$
\sup_t m_{\rm pole}(t)<-3.0484\times10^{-32},\qquad
\sup_t m_{\rm pole+details}(t)<-2.9782\times10^{-32}.
$$

For the35-channel result, the saved384-bit dyadic trial has
$-4.047\times10^{-32}<m(t_0)<-3.5939\times10^{-32}$,
paid gradient norm below $1.962\times10^{-17}$ and optimization
correction below $6.158\times10^{-33}$. Both results concern the
same complete retained direction and768 selected modes, under the
inherited actual common-form/domain and sound comparison/enclosure
premises. Neither comparison certifies this direction for any trial
on that support. This excludes neither arbitrary finer channels nor
the actual Schur form, other supports, full c27/cofinal positivity
or RH. No actual negative Weil test or new Lean result is supplied.

The new common contrast cache is
`gamma-rate-contrasts-depth1.json.zlib`. Reuse the saved exact trials
without another proposal solve:

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/gamma_support_optimization.py --finite docs/reports/weil-c27/gamma-support-6765-finite.json.zlib --keep-even-pole --trial docs/reports/weil-c27/gamma-support-pole-6765-optimum.json --out build/c27-even-pole-check.json
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/gamma_support_optimization.py --finite docs/reports/weil-c27/gamma-support-6765-finite.json.zlib --channel-refinement docs/reports/weil-c27/gamma-rate-contrasts-depth1.json.zlib --keep-even-pole --trial docs/reports/weil-c27/gamma-support-pole-contrasts-depth1-optimum.json --out build/c27-rate-details-check.json
```

The producer command below is for a genuinely missing common contrast
cache, not for consuming the archived result:

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/gamma_channel_resolution.py --out build/c27-rate-details.json.zlib
```

## Complete Gamma boundary on a common finite dual space

`full_gamma_dual.py` applies the existing inverse variational formula
and finite block elimination to `T=Lambda+Gamma_boundary+even_pole`,
with the same archived `Lambda=W^-1`, complete1281-coordinate retained
vector and selected support1281–2048. All Gamma rates are retained.
For a nonzero finite exact dyadic `w` outside that support, use the
common dual space `[E,w]`. Set

```
B0=E*TE, b=E*Tw, beta=w*Tw, a=w*r, z=w*Z, D=F-B0.
Q=v*Av-(F^-1 B_ret v)*F(F^-1 B_ret v), r=Cv-Z F^-1 B_ret v.
h=beta-z F^-1 z*+(z-b*)D^-1(z*-b).
sup_t U_[E,w](t)=Q-|a|^2/h >= sup_t m_T(t).
```

This is the supremum of a finite dual **upper bound**, not the full
comparison's exact optimum. Positivity of `D` follows from the
existing prime slack `D>=I/632`; the saved consumer additionally
checks768 directed positive pivots of `D-I/632` as finite numerical
consistency. This is not a new floor theorem. The two subsequent
FLINT solves each have one RHS. There is no new old-family optimizer.
The exact disjointness cancellation `z*-b=-E*C_p w` is available for
subsequent consumers.

The complete Gamma/pole kernel is separated from already cached
actual scalars by adding back every prime sine and cosine term and
subtracting the periodic Gamma diagonal. In particular its diagonal
keeps `-s_GP/(pi*n)`; removing just the cosine term is incorrect.
Digamma/trigamma full-rate formulas follow from
[DLMF5.7.6](https://dlmf.nist.gov/5.7.E6) and
[DLMF5.15.1](https://dlmf.nist.gov/5.15.E1). No finite-channel trace
allowance is applied to the full boundary. The streamed current-v
quadratic is new; the earlier support packet's `qv` belongs to an
older vector and is not relabeled. The original selected inverse
and transformed residual generators are reused without a new solve.

`full-gamma-dual-4096-finite.json.zlib` archives this new common action,
the exact384-bit dual vector on modes2049–4096 and the scalar inputs
needed to reconstruct its full boundary compression. Only the dual
trials are restricted; the inverse variational upper therefore needs
no debit for an omitted Fourier residual tail. Direct full-boundary
quadratics use independent diagonals and every cross term of `w`.

The directed result in `full-gamma-dual-4096-bound.json` is
**inconclusive**: its upper is approximately
`6.323379216251768e-31`, with actual current-direction
`Q≈2.971626487064166e-30`. A positive upper certifies no sign of the
full comparison or actual Schur form, and repairs no trial.

`full-gamma-dual-4096-diagnostic.json` uses the same exact `v,w,t0`
but the already studied35-channel/pole comparison. A small Woodbury
system gives its restricted fixed-trial upper
`4.847374044137042e-31>0`; the archived actual margin at that same
trial is below `-3.5939e-32`. Thus this dual restriction loses more
than `5.2067e-31` even before passing to all Gamma rates. It cannot
exclude that already excluded family. This identifies a limitation
of the chosen dual space; it does not determine full-Gamma
saturation or promise success from a larger cutoff.

The one production took14.23 seconds, the directed bound consumer
49.71 seconds and the diagnostic4.93 seconds at1024 working bits.
These are measured applications on the recorded environment, not
complexity or performance guarantees. The production reads the old
scientific scalar/inverse/residual packets and does not regenerate
any old target, Gram, scalar table or far moment. No new general
method, Lean verification, actual/cofinal positivity or RH result is
claimed. All numerical implications inherit the actual common-form,
domain and sound-comparison/enclosure premises of the earlier packet.

Consume the caches with fresh output paths:

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_dual.py --finite docs/reports/weil-c27/full-gamma-dual-4096-finite.json.zlib --out build/c27-full-gamma-dual-check.json
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_dual.py --diagnose --finite docs/reports/weil-c27/full-gamma-dual-4096-finite.json.zlib --out build/c27-full-gamma-dual-diagnostic-check.json
```

Only when these new common dual actions are absent, produce them with
`--produce --out <fresh-packet.json.zlib>`. Increasing working precision
in a cached consumer repeats no scientific production.

## Extend the same dual prefix through the complete archived middle

`full_gamma_dual_middle.py` produces only the new common actions for
a vector on2049–32768. Its old2049–4096 prefix is exactly unchanged,
including384-bit coefficients and original power-of-two normalization.
The current-direction Q and old prefix actions/scalar components are
reused. New residual and support-proposal rows start at4097, using
256/384-term finite geometric expansions; actual direct row checks
cover4097,16384,32768. Existing Cauchy projections give the new z and
u=z*-b=-E*C_p w, including their complete paid remainders. Only missing
Gamma/pole scalar components are separated from the archived actual
scalar table. No old actual scalar or far-moment supplier is rerun.

For the full-rate quadratic, set k=2pi/d, Phi_n=-omega_n*s_GP_n and
choose any common centering constant. On the finite supplied nodes,
write e_n=Phi_n-(pi/4)omega_n-constant. The linear term gives the signed
Hankel kernel J_mn=(-1)^(m+n)/(2(m+n)), including J_nn=1/(4n).
Finite scalar maxima give

```
eta=8 max|e_n|/[d*k^2*(2L+1)], L=2049;
delta=max|B_GP_nn-1/(4n)|;
|w*(B_GP-J)w|<=eta*(l1(w)^2-l2(w))+delta*l2(w),
l2(w)=sum|w_n|^2.
```

This is a finite triangle envelope from actual scalar enclosures,
with independent actual diagonals. The constant affects sharpness
only. No asymptotic formula, infinite-grid regularity, entrywise
Loewner claim or finite-channel trace allowance is assumed.

Existing FLINT polynomial multiplication of
`p(x)=sum(-1)^n*w_n*x^n` evaluates the real quadratic as
`sum_j coefficient(p*p,j)/(2j)`. Reuse the ACTUAL saved prefix beta,
add Lambda_tail and `Hankel_full-Hankel_prefix`, and pay the triangle
error on only the new offdiagonal pairs and tail diagonals. Both
prefix/tail cross terms are present. The exact saved prefix and its
parent hash guard this reuse; no changed normalization is hidden.
These are applications of existing polynomial convolution, finite
Cauchy sums and triangle inequalities, with no generic new algorithm.

`full-gamma-dual-32768-finite.json.zlib` has30720 exact dual coordinates
and31488 scalar rows1281–32768. It uses the same accepted finite-action
schema and unchanged `full_gamma_dual.py` consumers. Its beta allowance
is approximately2.27151e-5; wide ordered grid endpoints honestly retain
this allowance. The all-trial full-Gamma upper is enclosed between
2.88818e-31 and2.88820e-31, so it remains inconclusive. At the SAME
fixed t0 the35-channel control upper is approximately
-1.012916950382631e-32<0. This is a negative upper at that ONE fixed
trial; no all-t exclusion by this restricted bound is inferred.
Neither result gives a positive full-Gamma trial, an actual Schur sign,
full/cofinal positivity, Lean verification or RH resolution.

At1024 working bits the production took50.78s, bound consumer50.60s
and diagnostic35.74s. The first two figures include their complete
recorded actions, not a stand-alone multiplication benchmark or an
algorithmic complexity guarantee. The inherited actual common-form,
domain and sound-enclosure/comparison premises are unchanged.

Consume the new cache; do not rerun its producer:

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_dual.py --finite docs/reports/weil-c27/full-gamma-dual-32768-finite.json.zlib --out build/c27-full-middle-bound-check.json
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_dual.py --diagnose --finite docs/reports/weil-c27/full-gamma-dual-32768-finite.json.zlib --out build/c27-full-middle-control-check.json
```

If this new action is genuinely missing, produce it with
`full_gamma_dual_middle.py --out <fresh-packet.json.zlib>`.


## Complete finite dual optimization defect

`full_gamma_dual_defect.py` consumes the preceding accepted complete
middle vector and scalar bound. The complete $r_J$ and $M_Jw$ profiles,
and the two inverse vectors, were not saved by that bound. Only these
missing actions and their scalar defect are produced. Original $Q,v,w$,
$z,b$, all full-rate scalar rows and the original diagonal $W$ are reused.
No old channel optimizer, retained quadratic, scalar supplier or floor
check is rerun.

For the same coordinate embeddings $E,J$, $J=2049{:}32768$, let
$T=\Lambda+G$, $G$ the complete Gamma boundary plus even pole,
$F=E^*HE$, $D=F-E^*TE>0$, $Z_J=J^*HE$,
$B_J=E^*TJ$, $K_J=Z_J^*-B_J=-E^*C_{\rm p}J$ and

$$
M_J=T_J-Z_JF^{-1}Z_J^*+K_J^*D^{-1}K_J.
$$

The existing block factorization gives $M_J\succeq\Lambda_J>0$:
the block with upper row $(F,Z_J^*)$ and lower-right
$T_J+K_J^*D^{-1}K_J$ is $T_{[E,J]}$ plus
$[I;K_J^*D^{-1}]D[I,D^{-1}K_J]$. Minimize in the $E$ coordinate
and retain the $\Lambda_J$ part of the original diagonal lower bound.
For any exact dyadic $y$, square completion therefore yields

$$
L(y)=2\Re(y^*r_J)-y^*M_Jy,\quad e=r_J-M_Jy,
$$

$$
Q-L(y)-e^*W_Je\ \le\
\mathcal U_J:=Q-r_J^*M_J^{-1}r_J\ \le\ Q-L(y).
$$

This uses ordinary finite inverse variational, block and inverse-defect
estimates, not a new theorem. It applies to the entire current finite
space, with all inherited common-form/domain and sound-comparison premises.
A positive lower rules out a negative restricted upper in this finite
space. It does not give a positive lower for the actual comparison
supremum, a larger dual space or the actual Schur form. A large defect
upper by itself does not prove that the actual optimization gap is large.

The new adapter calls existing FLINT `arb_poly` multiplication for the
classical Toeplitz/Hankel squared-node action. It retains the independent
diagonals, parity, both mixed terms and complete cached Gamma/pole symbol.
No asymptotic or full weighted-trace estimate is used. The direct residual
uses the accepted GO/GP geometric action and absolute remainder. With
$p=F^{-1}z$ and $q=D^{-1}(z-b)$ from the same certified solver,

$$
M_Jw=\Lambda_Jw+J^*G(Jw-Eq)+J^*HE(q-p).
$$

The numerical packet `full-gamma-dual-32768-defect.json.zlib` saves all
$30720$ coordinates of each missing profile and both $768$-coordinate
inverse vectors. Taking $y=\vartheta w$ with the saved exact $512$-bit
dyadic scalar, the directed output gives

$$
2.48441\times10^{-31}<\mathcal U_J<2.88819\times10^{-31},
\qquad e^*W_Je<4.03773\times10^{-32}.
$$

Thus further direction optimization inside this same finite $J$ cannot
make its restricted upper negative. The larger-space and actual Schur
questions remain unresolved. The producer took $54.63$ seconds at $1024$
working bits and checked actual direct rows $2049,4096,32768$. This is an
application timing, not an algorithmic complexity claim. No Lean,
full/cofinal positivity or RH result is supplied.

To produce only the new profiles from the accepted scientific inputs:

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_dual_defect.py --out build/c27-full-middle-defect.json.zlib
```

Existing numerical kernels and APIs are described by
[Netlib's structured matrix multiplication](https://www.netlib.org/utk/people/JackDongarra/etemplates/node384.html),
[FLINT's ball polynomial API](https://flintlib.org/doc/arb_poly.html) and
[Python-FLINT 0.9.0's multiplication binding](https://github.com/flintlib/python-flint/blob/0.9.0/src/flint/types/arb_poly.pyx).

## Infinite residual payment for one positive primal direction

`full_gamma_primal_tail.py` consumes the accepted finite defect and
inverse profiles, matching current-direction far energy, original
constant far weight and the existing smooth symbol polynomial with
pointwise remainder decay. It runs no scalar or moment producer,
finite action, inverse solve, floor check or optimizer. The Fourier
cutoff remains $M_0=32768$.

The same inverse-residual square completion gives a positive sufficient
comparison after the complete residual has been paid. For finite $y$
outside $E$, put $Z=HE$, $Z_O=O^*Z$, $K=Z_O^*-E^*TO$ and

$$
p_y=F^{-1}Z_O^*y,\quad q_y=D^{-1}Ky,\quad
t_y=q_y-p_y,\quad g=-Eq_y+Oy,
$$

$$
My=T_Oy-Z_Op_y+K^*q_y,\quad
e=O(r_O-My),\quad L(y)=2\Re(y^*r_O)-y^*My.
$$

Then $Ft_y+Z^*g=0$, $r-Zt_y=Tg+e$, $E^*e=0$ and

$$
m_T(t_y)=Q-L(y)-e^*T^{-1}e
\ge Q-L(y)-e^*We.
$$

Only this finite-input action is needed, with bounded/coercive $T$,
square-summable $HE$ columns and the inherited actual domain.
No infinite minimax or bounded infinite $M$ assertion is used.
For the exact saved $y=\vartheta w$, the existing finite lower has
already paid $e_J^*W_Je_J$. Its actual inverse profiles define the
exact implicit correction $t_y=\vartheta(q-p)$; no midpoint rounding
or new solve occurs. The only additional debit is the disjoint
$m>M_0$ residual.

Let $A_0=M_0+1$, $k=2\pi/d$, $\gamma_0=d/\pi^2$ and

$$
c_0=\frac14-\frac{27^{-1/2}}{1-27^{-2}}+2(\cosh(d/2)-1),
\quad \rho_n=-kns_n^{\rm GP}-\pi kn/4-c_0.
$$

The accepted normalized polynomial has exact analytic $C_0=-\pi/4$.
Its pointwise remainder $\epsilon(m)=a(M_0/m)^{64}+b(M_0/m)^{65}$,
$a,b\ge0$, satisfies $m\epsilon(m)\le A_0\epsilon(A_0)$. Thus

$$
|\rho_m|\le R_\infty:=|-kM_0C_1-c_0|
+kA_0\left(\sum_{j=2}^{63}|C_j|(M_0/A_0)^j+\epsilon(A_0)\right).
$$

The $640$-bit coefficient grid and $832$-bit remainder grid are read
separately. The exported $C_0$ interval is never multiplied by
unbounded $m$. The smooth polynomial already contains the pole;
its prime sine terms are separate.

For $m>M_0\ge n>0$, the exact disjoint kernel is

$$
G_{mn}=(-1)^{m+n}\left[\frac1{2(m+n)}
+\gamma_0\frac{\rho_m-\rho_n}{m^2-n^2}\right].
$$

Set $s_n^2=((A_0-n)^{-2}+(A_0-n)^{-1})/(A_0+n)^2$.
The elementary bound $\sum_{m>M_0}(m^2-n^2)^{-2}\le s_n^2$
and weighted Minkowski give

$$
\mathcal B_G(x)=\sqrt{D_{\rm far}}\left[
\frac{\|x\|_1}{2\sqrt{M_0}}
+\gamma_0\sum_n|x_n|(R_\infty+|\rho_n|)s_n\right]
\ge\|G_{\rm far}x\|_W.
$$

It remains valid at $n=M_0$ and requires no separated geometric
ratio. For $x$ on $E$, reuse the old actual-symbol budget $B_*$:

$$
\mathcal B_H(x)=\sqrt{D_{\rm far}}\frac{2B_*}{\pi}
\frac{\|x\|_1}{\sqrt{M_0-2048}}
\ge\|H_{{\rm far},E}x\|_W.
$$

The effective action is
$G_{{\rm far},J}w+H_{{\rm far},E}(q-p)-G_{{\rm far},E}q$.
With the matching archived residual upper $E_r$,

$$
\|e_{\rm far}\|_W\le\sqrt{E_r}
+|\vartheta|\bigl(\mathcal B_G(w)+\mathcal B_H(q-p)
+\mathcal B_G(q)\bigr).
$$

`full-gamma-primal-tail.json` encloses the squared upper allowance
below $1.264968\times10^{-31}$. Subtracting it from the accepted
finite lower gives

$$
m_T(t_y)>1.2194448\times10^{-31}>0.
$$

Under the inherited actual form/domain, comparison and sound-cache
premises, the actual Schur quadratic in this one original saved
direction has the same strict lower bound, since $H^{-1}\preceq T^{-1}$.
This does not decide the full retained matrix or a cofinal family.
It neither proves RH nor supplies a new general theorem, algorithm
or Lean verification. A positive restricted upper alone had not
supplied this primal conclusion; the global residual payment is
essential.

The cache-only consumer took $0.36$ seconds at $1024$ working bits.
Reconsume the accepted data without regenerating any supplier:

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_primal_tail.py --out build/c27-full-gamma-primal-tail-check.json
```

## A common retained column from the complete Gamma direction

`full_gamma_common_column.py` supplies the missing joint projection

$$
b_{\rm fin}=R_J^*(W_J(Mw)_J-w)\in\mathbb R^{1281},
\qquad J=2049{:}32768.
$$

It consumes the accepted complete finite action, defect profiles,
primal tail and original residual generators. The existing disjoint
Cauchy product supplies all retained modes0–1280 in one adjoint
action, with512 terms and its absolute remainder. No scalar/moment
supplier, inverse solve, optimizer or old target production runs.
The saved $R_Jv_0$ is used only to check the new column's contraction;
it is not substituted for the whole retained-to-$J$ map.

The packet encloses the column on the640-bit grid at1024 working bits.
Its production took10.48 seconds. An independent1536-bit cached
consumption produced the same1281 column interval pairs. Direct
sums over all30720 exterior nodes at retained modes0,39,1280 lie in
those enclosures. The saved energy bounds are approximately

$$
E_M=14.1648613043077766231,\qquad
\delta_{\rm up}=298.7040165047722506653.
$$

Here $E_M$ is an upward square of the accepted unscaled effective
action tail norm. It is distinct from the tiny scaled residual
debit of the single-direction certificate. The denominator upper
is the finite original-$W$ energy plus $E_M$ minus the lower
endpoint of the same $h=w^*M_Jw$.

`full_gamma_common_target.py` applies the existing complete coupled
structured certificate to

$$
L_{\rm new}=L_0-\rho F_{\rm up}
+\frac{b_{\rm fin}b_{\rm fin}^*}{d_\rho},\qquad
d_\rho=\delta_{\rm up}+E_M/\rho>0.
$$

The inequality $L_{\rm new}\preceq S_N$ follows from the same
inverse-residual square completion and Young bound, with one common
linear trial for all retained coefficients. The full signed far
majorant $F_{\rm up}=P K P^*+\varepsilon I$ is reused. Its allowance
is read through the existing original-target loader on its512-bit
metadata grid; matrix entries keep their640-bit grid. The independent
diagonal is retained and updated separately. Signed $K$ is not
assumed positive, and no midpoint column outer product is used.
The new target includes all coupling to the1263-coordinate good
block; no old lift from a different target is reused.

At $\rho=1$, `full-gamma-common-rho1-pivots.json` certifies1263
positive pivots, then an inconclusive pivot, in54.82 seconds. The
existing `structured_residual_schur_finish.py` consumes this NEW
target and prefix to enclose its18-by18 complement in99.48 seconds.
All18 diagonal lower endpoints are positive; the full-spectrum
routine fails to isolate eigenvalues and retains none.

`small_schur_sylvester.py` instead calls FLINT's existing determinant
on leading principal submatrices of that saved enclosure. It stops
at order2: on remaining modes0,1 the determinant is strictly below
$-3.90918\times10^{-71}$. This rejects the $\rho=1$ sufficient lower
target, even though its reduced diagonals are positive. It supplies
no negative direction of actual $S_N$, no old-target inertia count,
and does not revoke the saved positive scalar direction. Each
determinant has a separate degree-dependent absolute export grid,
so a small determinant is not rounded onto the entry grid. The
first two minor tests took0.21 seconds; no spectral retry or
parameter scan is used.

Run only genuinely missing new data products; outputs must be fresh:

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_common_column.py --out build/c27-common-column.json.zlib
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_common_target.py --target-out build/c27-common-rho1-target.json.zlib --out build/c27-common-rho1-pivots.json
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/structured_residual_schur_finish.py --target docs/reports/weil-c27/full-gamma-common-rho1-target.json.zlib --pivots docs/reports/weil-c27/full-gamma-common-rho1-pivots.json --precision 1024 --out build/c27-common-rho1-schur.json
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/small_schur_sylvester.py --target docs/reports/weil-c27/full-gamma-common-rho1-target.json.zlib --pivots docs/reports/weil-c27/full-gamma-common-rho1-pivots.json --schur docs/reports/weil-c27/full-gamma-common-rho1-schur.json --out build/c27-common-rho1-sylvester.json
```

These are applications of existing numerical and variational methods,
with the inherited actual domain/comparison and sound-cache premises.
They supply no new generic algorithm/theorem or Lean verification.
Full retained positivity, cofinal positivity and RH remain unresolved.

## Reuse the old exact lift and determine its rank requirement

`full_gamma_lift_projection.py` consumes the OLD nineteen-channel
`joint-gamma-6765-rho32-v2-*` base, its1263-by18 saved inverse action,
the accepted common column and the original signed far factors. With
the old positive block denoted by (P_g), it uses the exact implicit
lift (V=[-P_g^{-1}B_g;I]) enclosed by that saved action. It supplies
only the missing products

$$
c=b_{\rm fin}^*V,\qquad \Beta_V=\Beta_0V,\qquad
G_V=(P_{\rm far}^*V)^*K(P_{\rm far}^*V)+\varepsilon V^*V.
$$

The signed112-feature core is retained, together with every cross
term in the independent allowance. Its exact rational coefficient
is the upper endpoint of the old base's explicitly640-bit transported
allowance, checked against the original512-bit metadata. Neither
the signed core alone nor an interval midpoint is assumed PSD.
This is the old base's complete lift, not a rounded saved direction
or the different lift of the failed common-rho1 target.

The1024-bit production took1.44 seconds. Independent1536-bit cached
consumption gave exactly the same3024 matrix interval pairs. No
target, inverse solve, channel values, special functions or moments
were regenerated.

`full_gamma_cross_screen.py` uses the saved channel inverse

$$
A_\rho^{-1},\qquad
A_\rho=I+A_0+(1+\rho^{-1})\delta_{\rm ch}I,\qquad \rho=1/32,
$$

from the last19-by19 block of the base kernel. This is not
(P_g^{-1}). Set (u=Mw), (F_\phi=u_J^*W_JU_J) and
(Z=A_\rho^{-1}\Beta_V). Since (A_0) includes all middle channel
rows, it majorizes the J-only Gram. Weighted Cauchy gives

$$
|(F_\phi Z)_j|\le\sqrt{e_J(Z^*A_0Z)_{jj}},\qquad
|(1+\rho)u_\infty^*W_\infty R_\infty V_j|
\le(1+\rho)\sqrt{E_M(G_V)_{jj}}.
$$

Subtract both radii from (|c_j|). All18 lower bounds are positive;
the smallest is above (1.8419798\times10^{-16}), at mode24. Thus
the full corrected common row has nonzero coordinates, under the
inherited common-object and sound-enclosure premises. This does
not supply a fixed positive matrix lower bound or pay coupling
to the good block. The screen took0.30 seconds; its90 interval pairs
match an independent1536-bit consumption exactly.

`small_schur_rank_floor.py` applies Jacobi's signature criterion to
the saved18-by18 complement of the SAME OLD nineteen-channel base.
All18 leading determinant intervals are strictly negative. Their
export grids are degree dependent, (640k+64), so the tiny higher
determinants are not rounded onto the entry grid. The pivot signs
are therefore one negative followed by17 positive. The complement
has inertia17 positive,1 negative,0 zero; with the1263 positive
prefix, the full sufficient base has1280 positive and1 negative
direction. This is not an inertia statement about actual (S_N).
A single positive update is permitted by this inertia count, but
its direction and paid strength still need verification. The18
determinant tests took0.25 seconds and were independently checked
on the saved small matrix, without an eigensolver or a large solve.

These consumers apply existing ball products, weighted Cauchy and
Jacobi/LDL signature; no new generic method or Lean verification
is claimed. Fresh-output usage:

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_lift_projection.py --out build/c27-lift-projection.json.zlib
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_cross_screen.py --out build/c27-cross-screen.json
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/small_schur_rank_floor.py --target docs/reports/weil-c27/joint-gamma-6765-rho32-v2-target.json.zlib --pivots docs/reports/weil-c27/joint-gamma-6765-rho32-v2-pivots.json --schur docs/reports/weil-c27/joint-gamma-6765-rho32-v2-schur.json --out build/c27-old-base-rank-floor.json
```

## Transport the scalar defect to the missing common cross row

`full_gamma_defect_cross.py` reuses information stronger than the
unscaled action energy. The accepted exact dyadic direction $v_0$
and exact nonzero $\vartheta$ have a certified finite defect

$$
e=R_Jv_0-\vartheta u_J,\qquad e^*W_Je\le E_e.
$$

Consequently the missing finite contraction satisfies

$$
F_\phi=(\Beta_0v_0)^*/\vartheta-e^*W_JU_J/\vartheta.
$$

Let $Z=A_\rho^{-1}\Beta_V$. The full common row has the enclosure

$$
\mathcal B_V\in
c-(\Beta_0v_0)^*Z/\vartheta+\mathcal E,
$$

$$
|\mathcal E_j|\le
\frac{\sqrt{E_e(Z^*A_0Z)_{jj}}}{|\vartheta|}
+(1+\rho)\sqrt{E_M(G_V)_{jj}}.
$$

The consumer reuses the first screen's finite radii with the
conservative factor $\sqrt{E_e/(\vartheta^2e_J)}$, taking the
numerator upward and $e_J$ downward. Its factor is less than
$0.109095699$. The far uncertainty remains paid. This obtains
common-row intervals from saved data without evaluating raw
channels, digamma functions, new moments or inverse actions.

Under the same actual form/domain and inverse-comparison premises,
the ordinary inverse-residual square completion can use

$$
\overline W=(W_J-W_JU_JA_\rho^{-1}U_J^*W_J)
\oplus(1+\rho)W_{\rm far}\succeq M^{-1},
$$

$$
\mathcal K=u^*\overline W u-h
=\delta+\rho\|u_{\rm far}\|_W^2
-F_\phi A_\rho^{-1}F_\phi^*\le
D_*:=\delta_{\rm up}+\rho E_M.
$$

Here $0<D_*<299.146668420532$.
The corresponding conditional monotone lower target is
$L_*+\mathcal B^*\mathcal B/D_*$. This uses the FULL common row;
it does not pair that denominator with the raw finite column.

All18 full-row intervals exclude zero. Mode0 is enclosed within
$(5.99690490660,9.65492947706)\times10^{-15}$. With the old core,
all18 quantities $(\Sigma_*)_{jj}+|(\mathcal B_V)_j|^2/D_*$
have positive lower endpoints; mode0's exceeds
$2.44652\times10^{-33}$. These are quadratic restrictions to the
old lift, not the Schur complement after updating its good block.
The production took0.41 seconds at1024 working bits.
Independent1536-bit cached consumption overlaps all128 interval
pairs;124 are endpoint-identical. All18 final row intervals,
diagonal lower bounds and metric budgets are endpoint-identical.
The four differences occur only in auxiliary center intervals.

For the full coupled test the scalar good-block metric
$\Omega=\mathcal B_GP_g^{-1}\mathcal B_G^*$ must still be enclosed
upward. An upper strictly below6.21431056 would preserve all18
updated diagonal lower bounds. It has not been supplied, and
positive diagonals would still leave the off-diagonal sign test.
The cached $P_g^{-1}B_g$ is not an inverse on this new row.
Coordinate interval boxes alone do not supply a fixed positive
matrix lower bound. Full retained/cofinal positivity and RH remain
unresolved; no general theorem, algorithmic novelty or Lean result
is claimed.

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_defect_cross.py --out build/c27-defect-cross.json
```

## Contract the shared row uncertainty against the positive complement

`full_gamma_correlated_sign.py` consumes the SAME OLD nineteen-channel
base, its accepted inertia, and the defect-transported full common row.
Write $b=\mathcal B_V$, and let
$\Gamma=\mathcal B_GP_g^{-1}\mathcal B_G^*\ge0$ be the exact missing
good-block metric. Standard bordered-matrix inertia gives, for $s>0$,

$$
\Sigma_*+b^*b/s\succ0
\quad\Longleftrightarrow\quad s+q<0,
\qquad q=b\Sigma_*^{-1}b^*.
$$

PSD corresponds to $s+q\le0$; exact equality leaves one zero eigenvalue.
For the exact updated lower target $s=D_*+\Gamma$. Replacing $\Gamma$
by a certified upper $\Omega$ gives a sufficient lower test; failure of
that weaker test need not exclude the exact update. These are existing
block-inertia/rank-one identities, not new general results.

Partition the saved core at its negative first pivot:

$$
\Sigma_*=\begin{pmatrix}a&r^*\\r&H\end{pmatrix},\qquad
P_+=H-rr^*/a\succ0,\qquad
T_+=\binom{-r^*/a}{I_{17}}.
$$

Its accepted one-negative inertia supplies $P_+\succ0$, without an
eigensolver or another determinant test. The scalar becomes

$$
q=|b_0|^2/a+(bT_+)P_+^{-1}(bT_+)^*.
$$

The consumer supplies the previously missing shared error Grams:
with $Z=A_\rho^{-1}\Beta_V$ and $\widehat b=c-(\Beta_0v_0)^*Z/\vartheta$,

$$
b=\widehat b+f+t,\qquad
f^*f\preceq C_f:=E_e\vartheta^{-2}Z^*A_0Z,\qquad
t^*t\preceq C_t:=(1+\rho)^2E_MG_V.
$$

Here $\rho=1/32$, $A_0$ includes all middle rows, and $G_V$ retains
the signed far core and the full allowance. Both inequalities concern
components of the same actual row, with the inherited common-form,
domain, inverse/comparison and sound-cache premises. No independent
attainment of error extrema is assumed.

One validated FLINT17-by17 solve with35 right-hand sides contracts
the center and both transported Grams. Define

$$
n_c^2=(\widehat bT_+)P_+^{-1}(\widehat bT_+)^*,\qquad
n_f^2=\operatorname{tr}(P_+^{-1}T_+^*C_fT_+),\qquad
n_t^2=\operatorname{tr}(P_+^{-1}T_+^*C_tT_+).
$$

Weighted trace domination and the triangle/reverse-triangle bounds give

$$
\max(0,n_c-n_f-n_t)^2
\le (bT_+)P_+^{-1}(bT_+)^*
\le(n_c+n_f+n_t)^2.
$$

All endpoints are directed, and the negative term is evaluated as
the interval $b_0^2/a$. The1024-bit consumption takes0.40 seconds.
Its present enclosure is too wide to determine the sign:

$$
-792<q<1.043\times10^{25},\qquad
-493<D_*+q<1.043\times10^{25}.
$$

These are outer enclosures, not attained extrema or evidence that the
actual positive energy is large. In particular, neither a positive
full-matrix metric budget nor failure of the exact lower update is
certified. Independent1536-bit consumption gives the same inconclusive
verdict: all648 Gram pairs and the four intervals for actual scalar
quantities overlap. Of all17 scalar exports,9 overlap and6 are identical;
the eight disjoint exports represent different computed bound values,
not conflicting enclosures of an actual quantity. The center-energy
width ratio is0.999999997454, so this comparison supplies no reason
to repeat the same calculation at higher working precision.

The packet retains both shared Grams and the17 scalar intervals.
Its center-only interval is already broad before either trace radius
is added. This calculation alone does not separate input uncertainty,
conditioning, lost joint cancellations and solve overestimation; no
condition number or spectral floor was measured. The exact-reference
certificate below supplies useful center and trace bounds from the
accepted packets. Their joint full-row sign and the good-block metric
remain separate obligations. Scalar noncancellation alone supplies
neither.
This applies existing ball solves, inertia and weighted trace bounds.
No large inverse, scientific producer, new generic theorem, Lean
verification, actual Schur-sign or RH conclusion is supplied.

Fresh-output usage:

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_correlated_sign.py --out build/c27-correlated-sign.json
```

## Certify the directional metric using an exact small reference

`full_gamma_directional_metric.py` applies the existing inverse-residual
identity and Neumann bound to the SAME actual positive complement
$P=P_+$ and $h=\widehat bT_+$. It reads the saved `finite_center`
fields, retaining their uncertainty; it does not substitute the full
row containing the finite/far errors.

Use exact rational midpoints of the mirrored old core to construct
$P_0,T_0$, and a rational reference $h_0$. FLINT computes the small
exact reference inverse $R=P_0^{-1}$ and $x=Rh_0^*$. Both identities
$RP_0=I$ and $P_0x=h_0^*$ are checked by exact rational products.
The midpoint is a reference, not an assumed actual or positive matrix.

Put $\Delta P=P-P_0$, $\delta h=h-h_0$ and $r=\delta h^*-\Delta Px$.
If $E\ge|\Delta P|$ entrywise and
$\eta\ge\||R|E\|_\infty$ is a certified upper strictly below one,
the existing residual identity gives

$$
hP^{-1}h^*
=h_0x+2\operatorname{Re}(\delta h\,x)-x^*\Delta Px+r^*P^{-1}r,
$$

$$
0\le r^*P^{-1}r
\le\frac{\|r\|_1\|Rr\|_\infty}{1-\eta}.
$$

Actual $P\succ0$ comes from the accepted old-core inertia. All source
uncertainty, including the negative pivot and $T_+$, remains enclosed.
The reference energy is computed exactly, signed perturbations are
kept, and the nonnegative residual allowance is added upward. No
cancellation-prone interval evaluation of $2hx-x^*Px$ is required.

The1024-bit certificate takes0.104 seconds and gives

$$
\eta<3.976\times10^{-28},\qquad
|2\operatorname{Re}(\delta h\,x)|<4.39\times10^{-79},
$$

$$
|x^*\Delta Px|<3.30\times10^{-28},\qquad
r^*P^{-1}r<1.309\times10^{-55}.
$$

Consequently

$$
391.85784253565945840673491942575
<hP^{-1}h^*<
391.85784253565945840673491942642.
$$

The directed interval width is below $6.584\times10^{-28}$.
Independent1536-bit consumption takes0.092 seconds: all324 exact
reference fields match, the actual correction/energy intervals
overlap, and both final widths satisfy that bound. These are exact
reference-arithmetic checks and numerical enclosures, not Lean proofs.
The reference profile and inverse have large norms, but no actual
condition number or intrinsic uncertainty floor is inferred.

`full_gamma_reference_trace.py` reuses the saved $R$ without an inverse
or solve. For $0\le t\le1$, $P_0+t\Delta P$ is nonsingular by the
same Neumann bound. Hermitian inertia is constant on this path;
actual $P\succ0$ therefore implies $P_0\succ0$. The eigenvalues of
$R\Delta P$, similar to a Hermitian matrix, have magnitude at most
$\eta$. Thus the standard Loewner comparison is

$$
\frac{R}{1+\eta}\preceq P^{-1}\preceq\frac{R}{1-\eta}.
$$

For each SAME actual PSD $C_j^+=T_+^*C_jT_+$, $j=f,t$,

$$
\frac{\operatorname{tr}(RC_j^+)}{1+\eta}
\le\operatorname{tr}(P^{-1}C_j^+)
\le\frac{\operatorname{tr}(RC_j^+)}{1-\eta}.
$$

The full trace includes all offdiagonal terms and the uncertainty in
$T_+$. The thin1024-bit contractions take0.030 seconds and yield
approximately $8.869400460878844$ and $6.026200142006464$; their
directed actual-trace interval widths are below $8.423\times10^{-27}$
and $5.210\times10^{-27}$ respectively.
Independent1536-bit thin consumption preserves the verdict and
overlaps the actual scalar intervals and final $q$ enclosure.
Of17 scalar exports,12 overlap and11 are identical; the five disjoint
exports are generated bound values.

Combining these with the certified center and saved full-row $b_0$
interval gives the simultaneous conservative enclosures

$$
-586<q<388,\qquad -287<D_*+q<688.
$$

These separate bounds leave the sign inconclusive. They do not
use every correlation between the negative coordinate and the
positive-subspace errors already recorded in $C_f,C_t$. The joint
signed estimate below uses these correlations to exclude the
$D_*$-based sufficient update. A sharper common gain charge and the
good-block inverse metric remain missing.
Neither interval endpoint is claimed attainable. All results retain
the inherited actual common-form/domain/inverse/comparison and
sound-cache premises. Existing exact rational, inverse-residual and
Loewner methods are reused; there is no new generic algorithm,
Lean-verification, actual full retained/cofinal positivity, Robin-tail
completion or RH claim.

Fresh-output usage:

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_directional_metric.py --out build/c27-directional-metric.json.zlib
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_reference_trace.py --out build/c27-reference-trace.json
```

## Use the negative/positive error correlations in the signed scalar

`full_gamma_joint_sign.py` reuses the same center profile, residual
energy ceiling $Q_c$, error traces $\tau_f,\tau_t$, and full shared
Grams. It evaluates only thin quadratic products. For the exact
center $\widehat b$, write

$$
q_c=\widehat b\Sigma_*^{-1}\widehat b^*
=|\widehat b_0|^2/a+hP^{-1}h^*,\qquad
g=\Sigma_*^{-1}\widehat b^*.
$$

Use the finite center coordinate here, not the full-row $b_0$ interval.
The accepted exact reference profile gives

$$
g=g_x+TP^{-1}r_c,\qquad
g_x=(\widehat b_0/a)e_0+Tx,\qquad
r_c=h^*-Px,\quad r_c^*P^{-1}r_c\le Q_c.
$$

For the SAME actual PSD $C_j$, trace domination bounds the correction
in its seminorm:

$$
\|C_j^{1/2}TP^{-1}r_c\|^2\le\tau_jQ_c,\qquad
L_j:=\sqrt{(g_x^*C_jg_x)_{\rm up}}+\sqrt{\tau_{j,\rm up}Q_{c,\rm up}}.
$$

Thus each corresponding linear error has magnitude at most $L_j$.
For $e=f+t$, the exact expansion is

$$
q(b)=q_c+2\operatorname{Re}(eg)+e\Sigma_*^{-1}e^*.
$$

The negative-pivot decomposition gives uniform quadratic-error bounds

$$
-\frac{(\sqrt{(C_f)_{00}}+\sqrt{(C_t)_{00}})^2}{-a}
\le e\Sigma_*^{-1}e^*
\le(\sqrt{\tau_f}+\sqrt{\tau_t})^2.
$$

All ceilings are directed, including the positive lower denominator
$(-a)_{\rm lo}$. The estimates retain the shared negative/positive
correlations in the linear contraction; they do not assume independent
error optimizers are simultaneously attained.

The1024-bit thin contraction takes0.021 seconds. It bounds the center
signed energy near $-128.17372946339$, the total linear error by
$62.805752720500$, and the quadratic error between
$-28.404847233170$ and $29.517335427025$. Consequently

$$
\boxed{-219.39<q<-35.85},\qquad
\boxed{79.76<D_*+q<263.30}.
$$

Independent1536-bit consumption takes0.017 seconds; all17 scalar
interval pairs and the final $q$ interval are endpoint-identical.

Under the inherited actual common-form/domain/inverse/comparison and
sound-cache premises, this excludes the PARTICULAR $D_*$-based
sufficient update for every $\Gamma\ge0$. It does not give a negative
direction for the actual Schur form or exclude stronger gain charges.

The negative upper for $q$ also supplies a concrete full-matrix budget.
If new certified positive $D_\sharp\ge\mathcal K$ and
$\Omega\ge\Gamma$ satisfy

$$
\boxed{D_\sharp+\Omega<35.85},
$$

the same block-inertia criterion makes that strengthened lower update
positive definite, including all coupled coordinates. These two bounds
are not supplied here. A known available refinement is

$$
\mathcal K\le D_*-F_\phi A_\rho^{-1}F_\phi^*.
$$

It is an inequality: $D_*$ is already a ceiling for
$\delta+\rho\|u_{\rm far}\|_W^2$. A quantitative common-source
channel-energy lower bound and the good-block metric remain missing.
The exact common charge must not be replaced by $D_*$ in an equality.
This is an application of existing quadratic expansion, trace and
weighted Cauchy methods, with no generic novelty or Lean verification.
Actual full retained/cofinal positivity, the Robin tail and RH remain
unresolved.

Fresh-output usage:

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_joint_sign.py --out build/c27-joint-sign.json
```

## A scalar-defect lower bound for the channel energy

`full_gamma_channel_charge.py` consumes the same old nineteen-channel
inverse, exact direction, scalar defect and joint signed packet. It
contracts the previously unused scalar channel energy; no inverse,
solve, scientific supplier or earlier center/trace/sign calculation runs.
Write

$$
F_0=(\Beta_0v_0)^*/\vartheta,\qquad
F_\phi=F_0-\eta,\qquad
\eta=e^*W_JU_J/\vartheta,\qquad k=E_e/\vartheta^2.
$$

Because $U_J^*W_JU_J\preceq A_0\prec A_\rho$, weighted Cauchy and
the operator norm imply

$$
\eta A_\rho^{-1}\eta^*\le k.
$$

For example the squared norm of
$W_J^{1/2}U_JA_\rho^{-1/2}$ is at most1. This bound uses an operator
norm, with no factor equal to the nineteen channels. Set

$$
E_0=F_0A_\rho^{-1}F_0^*,\qquad
L_F=(\max(0,\sqrt{E_0}-\sqrt{k}))^2,\qquad
U_F=(\sqrt{E_0}+\sqrt{k})^2.
$$

The triangle inequalities in the $A_\rho^{-1}$ norm give

$$
L_F\le F_\phi A_\rho^{-1}F_\phi^*\le U_F,\qquad
0\le\mathcal K\le D_*-L_F.
$$

The1024-bit computation takes0.413 seconds, including packet loading.
It gives $E_0$ near106.3370533152211, $\sqrt{k}<4.023460$,
$L_F>39.545559$ and $U_F<205.505007$. In particular the positive
certified ceiling

$$
D_\sharp=(D_*)_{\rm up}-(L_F)_{\rm low}<259.601110
$$

improves the coarse charge. Combining it with the already certified
$q$ upper gives a sufficient good-metric budget near-223.750468;
this separate estimate supplies no positive allowance for $\Gamma$.
That failure is inconclusive for $\mathcal K+\Gamma+q$: the same
$\eta$ occurs in both $F_\phi$ and the finite error of $b$, so a joint
estimate can retain cancellations lost by the two separate bounds.
The actual good-block inverse metric also remains missing.

The packet retains the new $F_0^*$ and $A_\rho^{-1}F_0^*$ columns for
subsequent common-source contractions. These are applications of
existing weighted Cauchy and norm inequalities, not new generic
methods. All conclusions inherit the same actual common-form/domain,
inverse/comparison and sound-cache premises. No Lean verification,
actual full retained/cofinal positivity, Robin tail or RH result is
claimed.

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_channel_charge.py --out build/c27-channel-charge.json
```

## A zero-error point obstructs the current error-envelope relaxation

`full_gamma_charge_center_obstruction.py` reads only accepted scalar
packets, using exact integer endpoint arithmetic. It runs no inverse,
contraction, scientific supplier or prior consumer. Its question is
whether the current *relaxed information* can force
$\mathcal K+\Gamma+q<0$; it does not seek an actual-source counterexample.

Keep the same cached $F_0,A_\rho,\widehat b,\Sigma_*$ and write

$$
b=\widehat b+\eta Z+t,\qquad F_\phi=F_0-\eta,
$$

$$
\eta^*\eta\preceq kA_0,\qquad t^*t\preceq C_t,\qquad\Gamma\ge0.
$$

These upper envelopes alone admit $\eta=t=0$. Put
$E_0=F_0A_\rho^{-1}F_0^*$ and
$q_c=\widehat b\Sigma_*^{-1}\widehat b^*$. At that relaxed point,

$$
D_*-E_0>192.809615,\qquad D_*-E_0+q_c>64.635885.
$$

Thus no uniform negative upper on this stated relaxation is possible
for the $D_*-F_\phi A_\rho^{-1}F_\phi^*$ envelope, even with $\Gamma=0$.
The stronger check also retains the same-source charge identity. Let
$J=\{2049,\ldots,32768\}$, $m_J=\|u_J\|_W^2$ and
$s=\|u_{\rm far}\|_W^2\ge0$. Then

$$
\mathcal K=m_J-h+(1+\rho)s-F_\phi A_\rho^{-1}F_\phi^*.
$$

The saved finite action energy and exact dual-denominator enclosure give

$$
m_J-h>284.539155,\qquad
m_J-h-E_0>178.202101,
$$

$$
m_J-h-E_0+q_c>50.028372.
$$

Consequently the relaxed point has $\mathcal K+\Gamma+q>50.028372$
for every $s\ge0$ and $\Gamma\ge0$. Even the ideal relaxed choice $s=0$
does not remove it. Sharpening only the far-energy upper or solving an
optimization problem over these unchanged upper envelopes cannot
supply the desired uniform negative bound.

The actual common source is fixed; $\eta,t,s,\Gamma$ also obey
relationships discarded by this relaxation. No assertion is made
that its zero-error point is realized by that source. The result
identifies missing information: use actual channel/cross actions,
further source constraints or a different lower construction. It
neither excludes the actual common gain nor proves actual Schur
negativity, full retained/cofinal positivity, a Robin tail or RH.
This is existing feasibility reasoning applied to cached numerical
scalars under the inherited common-form/domain, inverse/comparison
and sound-cache premises, without generic novelty or Lean verification.

```sh
python3 docs/reports/weil-c27/full_gamma_charge_center_obstruction.py --out build/c27-charge-center-obstruction.json
```

## The actual finite channel projection of the common source

`full_gamma_actual_channel_projection.py` supplies the missing
$F_\phi=u_J^*W_JU_J$ from the accepted cached $u_J=Mw$ and original
positive weights. The original nineteen channel shapes are evaluated
only for this new functional, one row at a time; their full values were
not archived. The effective action, residual, channel Gram and map,
exact lift and inverse actions are reused. No pole or contrast
rescaling is introduced. The channel definition source hash is checked
against the accepted boundary packet.

For the original bin $[a,b]$ and $m\in J$, the reused formula is

$$
U_{m,[a,b]}=(-1)^m
\sqrt{\frac{1-e^{-(4a+1)d/2}}{(b-a)d}}
\operatorname{Re}\left[
\psi_0\left(b+\tfrac14+i\pi m/d\right)
-\psi_0\left(a+\tfrac14+i\pi m/d\right)\right],
\qquad d=\log27,
$$

where $\psi_0$ denotes digamma. The new nineteen sums retain all
$m=2049,\ldots,32768$. With the cached $F_0$, they also supply
$\eta=F_0-F_\phi$ and the actual finite common row

$$
b_J=c-F_\phi A_\rho^{-1}\Beta_V.
$$

All eighteen coordinates overlap the separate transported expression
$\widehat b+\eta Z$. The1024-bit projection takes15.924 seconds and
gives channel energy near115.345991300848898, projected defect energy
near0.211760521889664071, and a positive common-charge ceiling below
183.800678. The earlier operator-only defect ceiling was above16.188229;
this projection supplies the actual finite channel information instead
of assigning that entire error envelope to every direction.

The packet preserves $F_\phi$, $\eta$, $A_\rho^{-1}F_\phi^*$ and $b_J$.
The packet supplies the finite channel data. The joint lower consumer
below controls the far terms through a common-source bound and uses
only nonnegativity of the good-block metric.
All conclusions inherit the actual common-form/domain,
inverse/comparison, channel-definition and sound-cache premises. This
is a new functional of the fixed source using existing channel formulas,
not a new generic method, Lean proof or retained/cofinal/RH result.

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_actual_channel_projection.py --out build/c27-actual-channel-projection.json.zlib
```

## The actual source excludes the fixed old-base rank-one lower update

`full_gamma_actual_joint_lower.py` uses the actual finite row $b_J$
above, the saved exact seventeen-dimensional reference inverse and
the accepted old far trace. Its new exact profile applies that inverse
to the new finite row; no inverse, solve, prior consumer or scientific
producer runs. The positive-complement residual certificate gives

$$
q_J=b_J\Sigma_*^{-1}b_J^*,\qquad
-120.742435<q_J<-120.742434.
$$

The new residual-energy ceiling is below $1.316\cdot10^{-55}$.
This row differs from the earlier scalar-defect center, so its signed
energy and gradient are newly contracted rather than copied from
the older center certificate.

Let $\kappa=1+\rho=33/32$, $s=\|u_{\rm far}\|_W^2\ge0$,
$E_\phi=F_\phi A_\rho^{-1}F_\phi^*$ and
$t=\kappa u_{\rm far}^*W_{\rm far}R_{\rm far}V$.
The same actual common source satisfies

$$
\mathcal K=m_J-h+\kappa s-E_\phi,\qquad b=b_J+t.
$$

Keep the full signed far Gram and its entire allowance, denoted $G$.
It majorizes $V^*R_{\rm far}^*W_{\rm far}R_{\rm far}V$.
For $g=\Sigma_*^{-1}b_J^*$, weighted Cauchy gives

$$
|t_0|^2\le\kappa^2sG_{00},\qquad
|tg|^2\le\kappa^2s\,g^*Gg.
$$

Using the accepted decomposition
$\Sigma_*^{-1}=a^{-1}e_0e_0^*+TP_+^{-1}T^*$, $a<0$, and
discarding only the nonnegative positive-complement part of the
quadratic far error yields

$$
\mathcal K+\Gamma+q
\ge m_J-h-E_\phi+q_J+A s-2B\sqrt{s},
$$

$$
A=\kappa-\frac{\kappa^2G_{00}}{-a},\qquad
B=\kappa\sqrt{g^*Gg},\qquad \Gamma\ge0.
$$

Here $A>0.7772979$. Completing the square in $\sqrt{s}$ gives
$As-2B\sqrt{s}\ge-B^2/A$ for every $s\ge0$.

To pay the new gradient without re-solving, use its new exact
reference profile $x_J$, residual $r_J=h_J^*-P_+x_J$ and
$g_J=(b_{J,0}/a)e_0+Tx_J$. The accepted far trace bounds

$$
\|C_t^{1/2}(g-g_J)\|\le\sqrt{\tau_tQ_J}.
$$

For the same fixed cap $E_M^0$ used in
$C_t=\kappa^2E_M^0G$, the saved positive lower $E_{M,\rm low}$
satisfies $E_M^0\ge E_{M,\rm low}$. Therefore

$$
B^2\le
\frac{(\sqrt{g_J^*C_tg_J}+\sqrt{\tau_tQ_J})^2}
{E_{M,\rm low}}.
$$

The1024-bit new contractions take0.108 seconds and give

$$
m_J-h-E_\phi+q_J>48.450729,\qquad B^2/A<10.153002,
$$

$$
\boxed{\mathcal K+\Gamma+q>38.297727.}
$$

Under the inherited actual common-form/domain, inverse/comparison
and sound-cache premises, this lower bound concerns the **actual
fixed common source**, not the earlier zero-error relaxation.
For every positive $D\ge\mathcal K$, it implies
$D+\Gamma+q>0$. The accepted one-negative-core criterion consequently
excludes PSD of the sufficient update
$L_*+\mathcal B^*\mathcal B/D$ for this fixed OLD base $L_*$ and
fixed saved dyadic dual $w$. The good-block inverse metric needs no
new computation for this exclusion.

The stronger base $Q-R^*\overline W R$ is not tested by that calculation
alone. Other dual directions and other lower constructions remain
outside this exclusion. A matrix
dominating this non-PSD lower bound can still be positive. This is
not actual $S_N$ negativity or an RH counterexample. Actual full
retained/cofinal positivity, the Robin tail and RH remain unresolved.
The methods are existing inverse-residual certification, weighted
Cauchy and square completion; no generic novelty or Lean proof is
claimed.

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_actual_joint_lower.py --out build/c27-actual-joint-lower.json.zlib
```

## The far-enclosure gap cannot rescue the same fixed dual

`full_gamma_baseline_gap.py` tests the specific stronger base
$\overline L=Q-R^*\overline W R$ with exactly the old nineteen
channels, their $A_\rho^{-1}$, $\rho=1/32$ and the saved fixed dual.
It reuses the accepted finite charge floor and exact seventeen-entry
trial. Only the changed small-matrix variational and far-gradient
contractions are new. No inverse, solve, old signed energy, channel
projection, action or trace supplier runs.

The exact baseline correspondence must distinguish the original
padding from its transported upper endpoint. Let $\widehat F$ be
the signed far approximation, $F_{\rm act}$ the actual far residual
Gram, $\epsilon_0$ the original complete padding, $\epsilon_1$ the
padding used by the later $\rho$ update, and $\epsilon$ the exact
rational upper in the accepted lift packet. Under the inherited
two-sided far enclosure,

$$
\|\widehat F-F_{\rm act}\|\le\epsilon_0\le\epsilon_1\le\epsilon.
$$

With $\kappa=1+\rho$ and
$C_J=Q-R_J^*W_JR_J+\Beta^*A_\rho^{-1}\Beta$, the source builders give

$$
L_*=C_J-\kappa\widehat F-(\epsilon_0+\rho\epsilon_1)I,
\qquad \overline L=C_J-\kappa F_{\rm act}.
$$

Consequently

$$
0\preceq\overline L-L_*\preceq2\kappa\epsilon I.
$$

There is no identification of $\epsilon_0$ with the transported
$\epsilon$. The original allowance metadata uses its explicit
512-bit grid; the reused upper has its separately labelled640-bit
grid. This comparison keeps the same $\overline W$ and does not
replace it by the different ideal $W_{\rm G}$ majorant.

For the accepted exact lift $V$, set $H_V=V^*V$ and

$$
\Sigma_{\max}=\Sigma_*+\delta H_V,\qquad
\delta=2\kappa\epsilon,
\qquad V^*\overline LV\preceq\Sigma_{\max}.
$$

The PSD addition preserves at least seventeen positive directions.
The new first pivot $a_{\max}$ is strictly negative, so this
eighteen-dimensional upper compression still has exactly one
negative direction and no zero direction. Its negative-pivot
positive complement $P_{\max}$ and lift $T_{\max}$ therefore give
the usual inverse decomposition. Reuse the saved exact trial $x_J$
from `full-gamma-actual-joint-lower.json.zlib`, without asserting
it is an inverse action for this changed complement. Define

$$
q_{J,\rm trial}
=\frac{|b_{J,0}|^2}{a_{\max}}
+2\Re(b_JT_{\max}x_J)-x_J^*P_{\max}x_J,
\qquad
g_x=\frac{\overline{b_{J,0}}}{a_{\max}}e_0+T_{\max}x_J.
$$

For the same actual $b=b_J+t$ and
$q_{\max}=b\Sigma_{\max}^{-1}b^*$, variational inversion gives

$$
q_{\max}\ge q_{J,\rm trial}+2\Re(tg_x)+|t_0|^2/a_{\max}.
$$

Use the full signed far Gram plus its entire allowance, denoted $G$,
with the same $t=\kappa u_{\rm far}^*W_{\rm far}R_{\rm far}V$
and $s=\|u_{\rm far}\|_W^2$. Weighted Cauchy and square completion
yield, if $A_{\max}>0$,

$$
\mathcal K+q_{\max}
\ge m_J-h-E_\phi+q_{J,\rm trial}-B_{\max}^2/A_{\max},
$$

$$
A_{\max}=\kappa-\kappa^2G_{00}/(-a_{\max}),\qquad
B_{\max}^2=\kappa^2g_x^*Gg_x.
$$

The1024-bit contractions take0.467 seconds. Directed endpoints give

$$
\delta<4.285\cdot10^{-135},\qquad
q_{J,\rm trial}>-122.995465,\qquad
B_{\max}^2/A_{\max}<10.153002,
$$

$$
\boxed{\mathcal K+q_{\max}>36.044697.}
$$

For every positive $D\ge\mathcal K$, strict positivity of
$D+q_{\max}$ excludes PSD of
$\Sigma_{\max}+b^*b/D$ by the existing one-negative-core criterion.
This matrix dominates the compression of
$\overline L+\mathcal B^*\mathcal B/D$, so that full stronger-base
update is also not PSD. A good-block inverse metric is unnecessary
for this exclusion.

The result is a conditional numerical application under the inherited
actual common-form/domain, inverse/comparison and sound-cache premises.
It excludes only this $\overline W$, saved fixed dual and admissible
positive denominators. Other majorants, new or combined directions
and other constructions remain outside its scope. Actual $S_N$ can
dominate the rejected lower bound and be positive. No actual $S_N$
negativity, full retained/cofinal positivity, Robin-tail, RH, generic
novelty or Lean proof is supplied.

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_baseline_gap.py --out build/c27-baseline-gap.json.zlib
```

## One new coupled two-source space

`full_gamma_two_source.py` and `full-gamma-two-source.json.zlib` extend
the preceding fixed-source exclusion to one specified two-dimensional
space. The old exact source $w_0$ and its actual effective action are
read unchanged. The new exact512-bit dyadic source $w_1$ is proposed
from $W_JR_JVe_0$, normalized by one power of two. The proposal is
rounded once and then treated as an arbitrary exact finite trial,
rather than the exact residual or an optimizer. Its integer-list
SHA256 is `cc753399b62b568f700b6fbcac4bd4816dbedf37d70d421561f8d7b89c15fec8`.
A nonzero exact two-coordinate minor certifies independence from $w_0$.
Both sources have support $J=2049{:}32768$.

Only $w_1$ needs new actual mixed actions, the corresponding two
certified768-dimensional solves, $u_1=Mw_1$, a1281-coordinate common
row and the original19-channel functional. The existing separated
Cauchy projection and FLINT Toeplitz/Hankel adapter supply these
actions. No old mixed solve, scalar table, old channel functional,
moment, residual generator or floor is regenerated. The new common
row is projected through every column of the saved exact18-column
lift. Direct actual rows at1281,2048,32768 and the two independently
computed selfadjoint crosses agree with their enclosures.

Write $U=[u_0,u_1]$, $W_{\rm src}=[w_0,w_1]$ and
$\Phi=U_J^*W_JU_{{\rm ch},J}$. With the same original19-channel
inverse and map $\Beta_V$, retain the full joint matrices

$$
K_J=U_J^*W_JU_J-W_{{\rm src},J}^*U_J
-\Phi A_\rho^{-1}\Phi^*,\qquad
B=(U_J^*W_J-W_{{\rm src},J}^*)R_JV
-\Phi A_\rho^{-1}\Beta_V.
$$

These include every source cross term. Use the preceding
$\Sigma_{\max}$, its negative first pivot $a$, positive complement
$P$ and elimination map $T$. Let $H=BT$, $B_0=Be_0$ and choose
one exact17-by2 rational trial $X$ with the saved reference inverse.
It is not an inverse action for the changed $P$. Put

$$
Q_X=B_0B_0^*/a+HX+X^*H^*-X^*PX,\qquad
g=e_0B_0^*/a+TX,
$$

$$
A=\kappa-\kappa^2G_{00}/(-a)>0,\qquad
C=K_J+Q_X-\frac{\kappa^2}{A}g^*Gg.
$$

Here $G$ is the same full signed far Gram plus its entire allowance.
For each coefficient vector $c$, the actual far row and far energy
belong to the common source $W_{\rm src}c$. The existing variational
inverse bound, weighted Cauchy and square completion therefore give
the uniform matrix inequality

$$
K_{\rm actual}+B_{\rm actual}\Sigma_{\max}^{-1}B_{\rm actual}^*
\succeq C.
$$

The far covariance is paid by this same-source inequality; it is not
set to zero or replaced by independently optimized scalar gains.
No numerical far profile of $u_1$ is required. The common domain
must contain the entire source span and both actions must have finite
weighted far energy; these remain explicit operator/model premises.

The1024-bit application takes77.55 seconds and gives, approximately,

$$
C=\begin{pmatrix}
36.0446979705376894&20.9502754574294758\\
20.9502754574294758&12.1909310344309248
\end{pmatrix}.
$$

The archived directed intervals, rather than this rounded display,
certify its first pivot $>36.044697$ and second Schur pivot
$>0.0139933305$. Thus $C\succ0$. For every positive-definite2-by2
denominator $D\succeq K_{\rm actual}$, standard bordered-block inertia
leaves one negative direction in
$\Sigma_{\max}+B_{\rm actual}^*D^{-1}B_{\rm actual}$.
Consequently the specified simultaneous two-source stronger-base
update is not PSD. The uniform inequality also excludes every
one-source combination from this span with its own admissible
positive denominator.

This conditional exclusion uses the same $\overline W$, support,
nineteen channels, $\rho=1/32$ and two-sided far-error comparison.
It excludes neither other dual spaces/majorants nor actual $S_N$
positivity. Full retained/cofinal positivity, Robin tail and RH remain
unresolved. This is an application with new actual-source data;
it supplies no new general theorem, algorithm or Lean verification.

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_two_source.py --out build/c27-two-source.json.zlib
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_two_source.py --precision 1536 --reuse-source docs/reports/weil-c27/full-gamma-two-source.json.zlib --out build/c27-two-source-independent.json.zlib
```

## Calibrate a source quadratic before selecting a new direction

`full_gamma_source_calibration.py` consumes only the accepted two-source
row/charge, eighteen-coordinate compression, far Gram and exact reference.
It tests a necessary small-matrix interface before constructing a whole-J
source eigensolver. No old source action, scalar/channel supplier, inverse
solve or optimizer runs.

Keep the preceding $a,T,P,G,A,\kappa$, and let $R_{\rm ref}$ be the saved
exact symmetric seventeen-dimensional reference. The new **linear** trial
rule is $X=R_{\rm ref}T^*B^*$, with

$$
L=e_0e_0^*/a+TR_{\rm ref}T^*,\qquad
Q=e_0e_0^*/a+T(2R_{\rm ref}-R_{\rm ref}PR_{\rm ref})T^*
-\frac{\kappa^2}{A}L^*GL.
$$

For one rational midpoint model, exact FLINT rational arithmetic checks
that the factored trial expression equals $K_J+BQB^*$. This is a different
trial rule from the archived fixed rational trial; their certificates
are not identified. The calibration scale is approximately36.04 and the
preselected tolerance is $10^{-8}\max(1,\|C_{\rm ref}\|_{\max})$.

| Evaluation on this same rational reference model | Measured maximum absolute error |
|---|---:|
| Ordinary doubles using the rounded **exact** $Q$ entries | $5.1006\times10^{89}$ |
| Factored doubles computing $R_{\rm ref}T^*B^*$ | $1.4422\times10^{161}$ |
| Factored doubles using rounded exact reference trial entries | $3.6005\times10^{87}$ |

All three fail the gate. No eigensolver or new source selection follows.
These are measurements for the recorded NumPy environment, not a
universal instability theorem or actual matrix condition number.

On the inherited640-bit input grid with1024 working bits, naive interval
evaluation of the variable linear trial has maximum output width
$1.7499\times10^{25}$; expanding through $Q$ gives
$2.7353\times10^{79}$. Fixing the exact rational reference $X$ first gives
width $1.4388\times10^{-27}$. **The last computation uses a different,
fixed trial rule**: it shares the rational midpoint value but does not
certify the variable-trial source operator. These widths measure the
specified evaluation, rather than intrinsic uncertainty. The frozen
trial calculation is an implementation option for a later certificate,
not a guarantee that a useful new source can be selected.

The canonical packet is `full-gamma-source-calibration.json.zlib`. It
records source/input hashes, the exact identity check, all three errors,
interval outputs, runtime versions and the failed calibration verdict.
This numerical diagnostic changes no accepted exclusion and settles no
new source space, actual Weil sign, cofinal positivity, Robin tail or RH.
It supplies no new general theorem, algorithm or Lean verification.

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 --with numpy==2.3.3 --with threadpoolctl==3.6.0 python docs/reports/weil-c27/full_gamma_source_calibration.py --out build/c27-source-calibration.json.zlib
```

## The common finite Gram on the exact eighteen-column lift

`full_gamma_lift_finite_gram.py` supplies the finite residual Gram missing
from the source-sign reduction. Keep the exact implicit lift $V$, its
eighteen-coordinate order, original nineteen channels and positive
diagonal weights. Put $Y=R_JV$, with $J=2049{:}32768$, and reuse
$\Beta_V=U_J^*W_JY$ and $A_\rho^{-1}$ from the accepted packets.
The new common products are

$$
G_J=Y^*W_JY,\qquad
S_J=G_J-\Beta_V^*A_\rho^{-1}\Beta_V=Y^*W_cY.
$$

The existing certified disjoint Cauchy action supplies the eighteen
new profiles with a fixed512-term expansion and its absolute geometric
remainder. The full30720-row Gram retains every cross term. Six direct
generator dot products check endpoint/interior rows in columns0 and17;
overlap is an implementation check, with provenance and the remainder
argument still required. No selected inverse, old source action,
scalar/channel supplier or complete1281-dimensional target is rerun.

The1024-bit application takes88.26 seconds on the recorded arm64,
Python3.13.12/Python-FLINT0.9.0 environment. The maximum exported entry
width for $S_J$ is below $4.132\times10^{-178}$. This is the width of
this chosen enclosure, rather than intrinsic input uncertainty or a
matrix condition number. The saved implicit lift is not rounded to a
new rational lift.

With $\kappa=33/32$, the same core, complete paid far Gram and exact
transported allowance also give

$$
\Xi=\Sigma_*+2\kappa\epsilon V^*V+\kappa G_{\rm far},\qquad
\Xi+\Theta\preceq\Xi+S_J,
\quad\Theta=Y^*(W_c-M_J^{-1})Y.
$$

The packet encloses entries of this exact **upper** matrix. Its scalar
entry endpoints are not Loewner matrix endpoints. A positive upper
matrix does not prove positive definiteness of $\Xi+\Theta$; no
threshold sign is tested here. The accepted far applicability condition
is reused, and no explicit $\Xi^{-1}$ is computed.

`full-gamma-lift-finite-gram.json.zlib` contains the five18-by18 matrices,
source/helper/input hashes, fixed expansion order, six direct checks,
grids, runtime and enclosure-width diagnostics. The transient full
profiles are represented by the pinned generators and lift, without a
second large archived copy. The common inverse action and complete
inverse-residual Gram remain missing. Actual source-domain and far-energy
identification, sufficient good-block/far certification, retained/cofinal
positivity, Robin tail and RH remain open. This is a numerical application
of existing Cauchy and matrix identities, without generic novelty or
Lean verification.

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_lift_finite_gram.py --out build/c27-lift-finite-gram.json.zlib
```

## Test the cached upper before requesting a common inverse

`full_gamma_threshold_upper.py` consumes only the saved eighteen-by-eighteen
matrix $U=\Xi+S_J$. The accepted same-object comparison gives

$$
\Xi\preceq H=\Xi+\Theta\preceq U,\qquad 0\preceq\Theta\preceq S_J.
$$

Since $\Xi$ has seventeen positive eigenvalues and $S_J\succeq0$, $U$
has at least seventeen positive eigenvalues. One existing Arb determinant
can therefore distinguish a positive upper from an upper with a negative
direction. This reuses finite Hermitian inertia, without a new generic
determinant method or another Gram/profile calculation.

The one1024-bit application takes0.53 seconds and encloses

$$
\det U\simeq1.06742140729487513\times10^{-1359}>0.
$$

Thus this exact upper matrix is positive definite under the inherited
premises. **This leaves $H$ undecided**: positivity of an upper does not
establish positivity of the actual threshold or a negative lower-source
direction. The inexpensive exclusion selector fails to exclude; this
route still needs a sharper common inverse/residual bound or another
sufficient ordered estimate.

`full-gamma-threshold-upper.json.zlib` records the outward determinant
endpoints, sign, input/helper/source pins and inherited conditions. The
matrix-entry grid remains640 bits; the determinant uses the existing
degree-based serialization rule $18\cdot640+64=11584$ absolute bits,
chosen before execution. There is no precision/grid scan, inverse/solve,
old supplier replay, source search, Lean verification, actual retained or
cofinal positivity, Robin-tail or RH result.

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_threshold_upper.py --out build/c27-threshold-upper.json.zlib
```

## Sharpen the upper with the two archived exact sources

`full_gamma_projection_upper.py` reuses $Z=[w_0,w_1]$ and their saved
$M_J$ actions on $J=2049{:}32768$. With the same exact eighteen-column
$Y=R_JV$, it supplies

$$
A_Z=Z^*M_JZ,\qquad B_Z=Z^*Y,\qquad
K_Z=B_Z^*A_Z^{-1}B_Z\preceq Y^*M_J^{-1}Y,
$$

and therefore the ordered upper

$$
\Xi\preceq H=\Xi+\Theta\preceq
U_2=\Xi+S_J-K_Z.
$$

The source Gram reuses its archived first diagonal and both directed
cross enclosures. Only its unexported second diagonal is contracted from
the cached action. One existing two-column Cauchy projection supplies
the missing $Z^*R_J$, then multiplies the same implicit $V$. Its input
contains unweighted exact sources, with the leading768 rows zero; the
original residual coefficient and zero-mode normalization are preserved.
Six direct generator-column checks cover both sources and retained modes
0,39,1280. A native two-by-two solve with eighteen right sides forms all
cross terms of $K_Z$.

The fixed1024-bit application takes9.17 seconds and encloses

$$
\det U_2\simeq-2.2564304658181403\times10^{-1361}<0.
$$

A1536-bit execution of the same consumer takes10.55 seconds. All724
matrix-entry interval pairs have identical endpoints; the determinant
intervals overlap with nonidentical endpoints and remain strictly
negative. Both executions pass the six direct checks, with all40
guarded files unchanged. This is a reproduction of the same algorithm.

Under the inherited channel-aware inverse comparison and seventeen
positive directions of $\Xi$, this gives seventeen positive and one
negative direction of $U_2$, and hence of $H$. The accepted inertia
identity gives $n_-(C_{\rm sharp})=0$ for this whole finite-$J$,
fixed-eighteen-lift lower source comparison, including sources outside
$\operatorname{span}Z$. Thus its strict-negative-source criterion is
excluded without supplying the complete common inverse. This does not
exclude actual sufficient updates generally or settle semidefinite
endpoint updates.

`full-gamma-projection-upper.json.zlib` saves $A_Z,B_Z,A_Z^{-1}B_Z,K_Z,U_2$,
directed pivots, determinant endpoints, projection checks and source/input/
helper pins. Entry intervals enclose the same exact matrices; their
endpoints are not Loewner bounds. The512-term expansion,640-bit entry
grid and11584-bit absolute determinant grid are fixed. No forward lift
profiles, old source actions, scientific suppliers, F/D solves or full
Gram are regenerated. Actual common-domain/far-energy identification,
sufficient good-block/far certification, retained/cofinal positivity,
Robin tail and RH remain open. This is a numerical application of
existing variational projection and inertia, with no generic novelty
or Lean-verification claim.

```sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_projection_upper.py --out build/c27-projection-upper.json.zlib
```

## A necessary quadratic for one archived-source surplus shift

[The program and canonical packet](https://github.com/the-omega-institute/trureturing-experiments/tree/ca237495af332048e91aeb50d5be69480014105a/docs/reports/weil-c27)
consume the original exact384-bit $w_0$. This changes the effective
comparison using actual surplus; it does not reopen the unchanged
comparison's source search.

Keep the SAME actual exterior $H$, complete $T=\Lambda+B_{\rm GP}$,
selected support $E$ and common effective action $\mathcal M$. With
$O$ its complement, $D=E^*(H-T)E$ and $C=O^*(H-T)E$, the actual
selected exterior Schur form has the conditional decomposition

$$
\mathcal A=\mathcal M+\Pi,\qquad
\Pi=(H-T)_{OO}-CD^{-1}C^*\succeq0.
$$

Fix the SAME full $u=\mathcal Mw_0$, $h=w_0^*u>0$. If a directional
budget $uu^*\preceq p\Pi$ and $p\le h(1-t)/t$ hold for $0<t<1$, then

$$
\mathcal M_t=\mathcal M+\frac{t}{h(1-t)}uu^*\preceq\mathcal A,\qquad
\mathcal M_t((1-t)w_0)=u.
$$

Thus the full charge and row change by
$\mathcal K_t=\mathcal K+th$ and
$\mathcal B_t=\mathcal B+t\,w_0^*R_{\rm res}$. The actual $u_{\rm far}$
is unchanged, not multiplied by $1-t$. No such directional budget
is certified by this consumer.

For the SAME old nineteen-channel base, positive good block,
actual complement $\Sigma_*$ and exact lift $V$, the actual sufficient
score includes its new good-block charge $\Gamma_t\ge0$.
The inherited order $\Sigma_*\preceq\Sigma_{\max}$, together with
their same nonsingular inertia, gives
$\Sigma_{\max}^{-1}\preceq\Sigma_*^{-1}$. The existing common far
square completion, with its complete signed/paid $G$ and
$\Xi=\Sigma_{\max}+\kappa G$, then supplies the necessary LOWER gate

$$
\mathfrak F_{w_0}(t)
=\mathcal K_t+\Gamma_t+
(\mathcal B_tV)\Sigma_*^{-1}(\mathcal B_tV)^*
\ \ge\ \ell(t)
=k_J+th+(b_J+t r_V)\Xi^{-1}(b_J+t r_V)^*.
$$

This retains the actual common forms/domains/inverses, same
source/lift/base, finite unscaled weighted far energy and sound
enclosure premises. A negative $\ell$ would not certify a sufficient
update. A nonnegative minimum would stop only this fixed-source
strict criterion; zero thresholds alone do not exclude shifted PSD
endpoints.

The three accepted inputs already contain all required data:

| Packet | Fields reused |
|---|---|
| full-gamma-two-source.json.zlib | finite_charge_matrix[0][0] is $k_J$; finite_row_on_lift[0] is $b_J$ |
| full-gamma-projection-upper.json.zlib | matrices.source_gram[0][0] is $h$; matrices.source_on_lift[0] is $r_V$ |
| full-gamma-lift-finite-gram.json.zlib | matrices.far_shifted_core is $\Xi$ |

[full_gamma_surplus_shift_gate.py](full_gamma_surplus_shift_gate.py)
validates those exact packet pins and their common metadata.
One native18-by18 solve with two right sides supplies the common
contractions and the three coefficients of $\ell(t)$. Exact rational
minima of the lower- and upper-coefficient polynomials enclose
$\min_{[0,1]}\ell$. Here $t,t^2\ge0$ justify scalar coefficient
endpoints; no matrix endpoint ordering is used. The closed interval
bounds the polynomial, not an operation at $t=1$.

The1024-bit run takes0.0657 seconds and returns the enclosure,
displayed approximately,

$$
-4.03437\times10^{26}
\ \le\ \min_{0\le t\le1}\ell(t)
\ \le\ 4.60803\times10^{24}.
$$

Its verdict is **inconclusive**: neither a negative region nor
nonnegativity is certified. These are widths from this evaluation,
not measured condition numbers, intrinsic uncertainty or a proof
that the shift cannot work. No precision retry or second scientific
execution follows this inconclusive gate. The different-working-directory
and path-with-spaces execution has not been tested.

[full-gamma-surplus-shift-gate.json.zlib](full-gamma-surplus-shift-gate.json.zlib)
saves the coefficient intervals, 43 exported intervals in total,
exact rational minimum bounds, solve result and provenance.
No old source/action/Gamma/Gram/projection production, new source,
large inverse or parameter scan is performed. A useful correlated
or factored certificate for this same minimum remains missing.
Any later negative necessary region still needs a SAME-$u$
directional surplus budget and an actual joint sufficient upper
including $\Gamma_t$, $\Sigma_*$ and the full far response.
No actual retained/cofinal positivity, Robin-tail, RH, generic novelty
or Lean result is supplied.

~~~sh
uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_surplus_shift_gate.py --out build/c27-surplus-shift-gate.json.zlib
~~~

## An admissible surplus interval from the same cached action

[full_gamma_surplus_budget.py](full_gamma_surplus_budget.py) supplies
the directional budget required by the preceding shift, under the
inherited actual form/domain and sound-cache premises. It reuses
the original unscaled $u=\mathcal Mw_0$; no action or energy is rerun.

The complete identity $H-T=A_\Gamma-\Lambda-C_{\rm p}$ and the
existing global weighted-prime slack give $H-T\succeq I/632$.
Indeed $q(59/10-r_q)>1/400$, $q\le79/50$ and
$\Lambda\preceq A_\Gamma-(59/10)I$ apply on the common form domain.
Taking the infimum over the selected finite $E$ coordinate yields
$\Pi\succeq I/632$ on $O$; a bound only on $D$ would not give this.

The common-column packet already contains an enclosure of the SAME
full action's middle weighted energy $E_J$ and an upper bound for its
far weighted energy $E_F$.
The far weight is the original positive CONSTANT
$W_{\rm far}=D_{\rm far}I$. The pilot supplies $D_{\rm far}$ and
all middle weights. With $w_{\min,J}=\min_{2049\le m\le32768}W_m$,

$$
\|u\|_2^2\le U_u:=\frac{E_J^{\rm up}}{w_{\min,J}}
+\frac{E_F^{\rm up}}{D_{\rm far}},\qquad
uu^*\preceq U_uI\preceq p_{\rm up}\Pi,\qquad p_{\rm up}=632U_u.
$$

For the accepted directed $0<h_{\rm low}\le h$, define

$$
t_{\rm cap}=\frac{h_{\rm low}}{h_{\rm low}+p_{\rm up}}.
$$

Every fixed rational $0<t\le t_{\rm cap}<1$ then satisfies
$p_{\rm up}\le h(1-t)/t$. This admits the stated
$\mathcal M_t\preceq\mathcal A$ and
$\mathcal M_t((1-t)w_0)=u$, while retaining the original unscaled far
action. It is a sufficient interval, not a maximization of $t$.

| Accepted packet | Fields consumed |
|---|---|
| full-gamma-common-column.json.zlib | bounds.middle_action_energy and bounds.far_action_energy_upper |
| inverse-residual-pilot.json.zlib | middle_weight_numerators, middle_weight_bits and the exact constant far weight |
| full-gamma-dual-32768-finite.json.zlib | original exact384-bit source and common metadata |
| full-gamma-dual-32768-defect.json.zlib | same action/source correspondence |
| full-gamma-primal-tail.json | provenance of the already paid weighted far energy |
| full-gamma-projection-upper.json.zlib | matrices.source_gram[0][0] for $h$ and the same row0 source hash |

The middle array starts at1281, so the $J$ slice begins at index768
and has30720 entries. The program validates the exact packet pins,
metadata edges and original source hash. All new arithmetic uses
standard-library exact rational numbers; no precision parameter,
solve, source/action, Gram, projection, tail producer or parameter
scan is invoked. The output stores every scalar as a rational pair
and one interior dyadic parameter, and refuses to overwrite an
existing output. The new scalar conclusion is conditional on the
actual common-form identification and accepted energy enclosures;
those premises are not discharged by reading the archives.

The canonical [full-gamma-surplus-budget.json](full-gamma-surplus-budget.json)
records the scalar bounds and parameter. The exact rational outputs are displayed approximately below;
the JSON pairs, rather than these rounded values, define the budget.

| Quantity | Approximate value |
|---|---:|
| $U_u$ | $4541.41732972820$ |
| $p_{\rm up}$ | $2870175.75238822$ |
| $t_{\rm cap}$ | $0.000374611013884066$ |
| Eligible exact interior example | $2^{-13}=1/8192$ |

The sole exact scalar application takes0.1789 seconds. Its terminal
verdict is conditional-positive-same-action-budget-interval. No
second scientific execution or precision retry is performed.

The earlier necessary-quadratic packet remains inconclusive and is
not recomputed. Eligibility of a shift supplies no sign of $\ell$
and no actual joint sufficient upper including $\Gamma_t$, the SAME
actual complement and full far response. No actual retained/cofinal
positivity, Robin-tail, RH, generic novelty or Lean result follows.

```sh
python3 docs/reports/weil-c27/full_gamma_surplus_budget.py --out build/c27-surplus-budget.json
```

## A cached fixed-trial exclusion over the sufficient surplus interval

[full_gamma_small_shift_exclusion.py](full_gamma_small_shift_exclusion.py)
uses the stronger-base compression already defined above, with the
SAME nineteen channels, $\rho=1/32$, $\overline W$, original exact384-bit
$w_0$, lift $V$ and full unscaled $u=\mathcal Mw_0$. It tests the
entire saved sufficient interval $0\le\tau\le t_{\rm cap}$; no new
parameter, source or cutoff is selected.

Let $a=a_{\max}<0$, $T=T_{\max}$ and $x_J$ denote the accepted
negative-pivot decomposition and saved exact seventeen-vector for
$\Sigma_{\max}$. This is the upper compression of the stronger
base $\overline L$, not $\Xi$ or the actual old $\Sigma_*$. Put

$$
r=w_0^*R_JV,\qquad g_0=Tx_J+(\overline{b_{J,0}}/a)e_0.
$$

The same fixed trial's signed value at $b_J+\tau r$ differs from
its old value by

$$
2\tau\Re(rg_0)+\tau^2|r_0|^2/a.
$$

Its far gradient is $g_\tau=g_0+\tau(\overline{r_0}/a)e_0$.
Reuse the accepted square-completion coefficient $A_{\rm low}>0$,
old penalty upper $P_0$ and old lower $L_0$ verbatim. For the SAME
complete paid PSD far Gram $G$, choose a directed upper

$$
P_1\ge\kappa^2|r_0/a|^2G_{00}/A_{\rm low}.
$$

The $G$-metric triangle inequality bounds the new completed-square
penalty by $(\sqrt{P_0}+\tau\sqrt{P_1})^2$. Thus

$$
\mathcal K_\tau+q_{\max}(\mathcal B_\tau V)
\ge L_0+c_1\tau+c_2\tau^2,
$$

$$
c_1\le h_{\rm low}+2\Re(rg_0)-2\sqrt{P_0P_1},\qquad
c_2\le |r_0|^2/a-P_1.
$$

Here $\mathcal K_\tau=\mathcal K+\tau h$,
$\mathcal B_\tau=\mathcal B+\tau w_0^*R_{\rm res}$, and the original
far action is unscaled. The fixed-trial lower, actual common forms,
comparison domains, complete far enclosure, same-source identity
and sound cached intervals remain inherited conditions.

| Accepted packet | Reused fields |
|---|---|
| full-gamma-baseline-gap.json.zlib | shifted negative pivot and gap coefficient, $A_{\rm low}$, $P_0$, $L_0$ |
| full-gamma-actual-joint-lower.json.zlib | new_reference_profile, used as an exact fixed trial only |
| joint-gamma-6765-rho32-v2-schur.json | reduced_schur_intervals, first row only |
| full-gamma-lift-projection.json.zlib | lift_gram first row and complete far_gram_upper diagonal00 |
| full-gamma-actual-channel-projection.json.zlib | actual_finite_row_on_lift, coordinate0 |
| full-gamma-projection-upper.json.zlib | matrices.source_on_lift[0], the same source's $r$ |
| full-gamma-surplus-budget.json | exact h_lower and admissible_parameter_upper |

The program reconstructs only the needed first-row entries of
$\Sigma_{\max}$, evaluates the new mixed scalar $rg_0$ and shifted
far-penalty scalars, then exports downward coefficient endpoints.
It does not repeat the old trial quadratic, old gradient energy,
old penalty/lower, any source/action, inverse, solve, Gram, far
producer or reference profile. The existing exact rational
quadratic-minimum helper is reused after substituting
$\tau=t_{\rm cap}s$, $0\le s\le1$. All input bytes are hash-pinned
before and after consumption; output is exclusive and never overwrites
an accepted packet. A failed prerequisite or nonpositive lower minimum
stops without an automatic supplier, precision or parameter retry.

The sole1024-bit application takes0.0113 seconds. Its640-bit outward
enclosures and exact rational coefficients are in
[full-gamma-small-shift-exclusion.json](full-gamma-small-shift-exclusion.json).
The following displays are approximate; the JSON rational pairs define
the actual lower polynomial and its minimum.

| Quantity | Approximate value |
|---|---:|
| $L_0$ | $36.0446979705377$ |
| $rg_0$ | $-954.862843927571$ |
| $P_1$ | $8005.06261257558$ |
| $c_1$ | $-1404.30017560541$ |
| $c_2$ | $-32507.0069794674$ |
| Minimum of the lower polynomial | $35.5140698387622$ |
| Minimizing parameter | $t_{\rm cap}\approx0.000374611013884066$ |

Under the stated premises, this gives
$\mathcal K_\tau+q_{\max}(\mathcal B_\tau V)>35.514069$ throughout
the interval. For every eligible rational $0<\tau\le t_{\rm cap}$
and every positive $D\ge\mathcal K_\tau$, the one-negative-core
criterion excludes PSD of
$\Sigma_{\max}+(\mathcal B_\tau V)^*(\mathcal B_\tau V)/D$.
Since this matrix dominates the compression of
$\overline L+\mathcal B_\tau^*\mathcal B_\tau/D$, that specified
lower update is also excluded. No1263good-block inverse upper is
needed for this exclusion. The zero endpoint is the old construction.

This leaves the earlier inconclusive $\ell$ packet unchanged and does
not turn its lower relaxation into an actual signed-inverse upper.
The budget interval is only sufficient; larger truly eligible shifts,
other sources/bases/majorants and other constructions remain unexcluded.
The rejected lower update can be dominated by a positive actual Weil
operator. No actual Weil negativity, retained/cofinal positivity,
Robin-tail, RH, generic novelty or Lean result is supplied.

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_small_shift_exclusion.py --out build/c27-small-shift-exclusion.json
```

## The saved prime-potential step weight cannot meet the fixed1/48 budget

[full_gamma_prime_potential_budget.py](full_gamma_prime_potential_budget.py)
reuses the original weighted row and the complete source profile. Put
$v=59/10-r_q$ and keep the SAME exterior form, $E$, $O$, full unscaled
$u=\mathcal Mw_0$ and $h=w_0^*u>0$. Under the inherited common-domain
identifications, $H_{\rm ext}-T\succeq\operatorname{mult}(v)$.
For every admitted $f=e+y$, orthogonality and weighted Cauchy give

$$
|\langle u,y\rangle|^2\le p_v(H_{\rm ext}-T)[e+y],\qquad
p_v=\int_I|u|^2/v.
$$

Taking the infimum over the finite $E$ coordinate gives
$uu^*\preceq p_v\Pi$ for the SAME closed short. Multiplication by $v$
need not preserve either Fourier subspace. This uses standard weighted
Schur, Cauchy and form-short estimates; no new general theorem is claimed.

The existing `high-floor.json` supplies26disjoint outward cells on
$s=x/L\in[0,1]$. Its saved minimum $\mu_C$ lower-bounds $qv$,
so the reused cell floor is

$$
v_C=\frac{\mu_C}{1+(29/50)(\sup C)^2},\qquad
W|_{\pm LC}=v_C^{-1}\ge v^{-1},\qquad m=\max_Cv_C^{-1}.
$$

No old logarithm, prime, source/action, Gamma, inverse, solve, Gram
matrix or far producer is repeated. The cached source uses the original
unphased orthonormal cosine basis. FLINT polynomial convolution and
autocorrelation retain both sum and difference frequencies through65536.
All input coefficients remain interval balls. Rational-angle sine/cosine
anchors and the addition formula in blocks of64 integrate the square on
each saved cell at1024working bits, with640-bit outward exports. Exact
endpoint primitives at0and1 retain the Parseval normalization. The new
energy computation takes1.264seconds and has no parameter/cell/precision
scan. The complete unknown far and its cross term are paid by the weighted
triangle; no far phase or cancellation is assumed.

| Quantity | Approximate value |
|---|---:|
| $\min_Cv_C$ | $0.00174356661223271$ |
| $\|u_J\|^2$ | $1344.75921630745$ |
| $a_1=\langle u_J,Wu_J\rangle$ | $656214.527747263$ |
| $mU_F$ | $26846.3715529475$ |
| $(\sqrt{a_1^{\rm up}}+\sqrt{mU_F})^2$ | $948518.835933079$ |
| $47h$ | $50553.3119129287$ |

The forward upper does not certify the preselected $\tau=1/48$.
The exact data and all26cell energies are in
[full-gamma-prime-potential-budget.json](full-gamma-prime-potential-budget.json).
They also permit a stronger stop test without another Fourier calculation.

[full_gamma_prime_potential_obstruction.py](full_gamma_prime_potential_obstruction.py)
uses only exact rational postprocessing of these cells and the original
complete far cap. Let $a_-$ be the sum of nonnegative lower cell-energy
endpoints times $W_C$. Integer square roots certify

$$
\sqrt{a_-}\ge810,\qquad \sqrt{mU_F}\le164.
$$

The reverse weighted triangle therefore gives, for the SAME full source
and every far phase satisfying its complete cap,

$$
\boxed{\int_IW|u|^2\ge(810-164)^2=417316>47h.}
$$

The comparison uses the archived upper endpoint for $h$, not its lower.
The scalar-only application takes0.005seconds; its exact certificate is
[full-gamma-prime-potential-step-obstruction.json](full-gamma-prime-potential-step-obstruction.json).
Thus this exact saved26-cell step weight cannot supply the required
multiplication-weight budget at1/48. More accurate integration or a sharper
cross estimate with this same $W$ cannot remove that obstruction.

This lower bound is for $\int W|u|^2$, not for $p_v\le\int W|u|^2$.
It neither rejects the exact reciprocal potential, a different majorant,
the optimal $\Pi$ budget nor actual Weil positivity. Actual common forms,
source/far identity and sound saved intervals remain conditions. Eligibility
with a different estimate, the actual good-block sufficient score, cofinal
positivity, Robin tail and RH remain unresolved. No Lean or originality
claim is supplied.

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_prime_potential_budget.py --out build/c27-prime-potential-budget.json
PYTHONDONTWRITEBYTECODE=1 python3 docs/reports/weil-c27/full_gamma_prime_potential_obstruction.py --out build/c27-prime-potential-step-obstruction.json
```

The scalar consumer pins the published energy packet; its default input
directory is the report directory. Both outputs must be fresh paths.

## A refined reciprocal-potential budget at the fixed shift1/48

[full_gamma_prime_potential_refined.py](full_gamma_prime_potential_refined.py)
keeps the same full source, actual comparison/short and complete far cap.
It reuses the first25 saved cell energies and weights, the saved norm
$a_0$, and the total energy of the last cell. Only the last cell
$[l,1]$, $l=261826231978/274653072167$, is divided into128 equal rational
subcells, fixed before the sole integration run. On each new $[a,b]$,

$$
\widetilde W_C=
\frac{1+(29/50)b^2}{\min_{[a,b]}P},
\qquad 1/v\le\widetilde W\le W.
$$

Here $P\le qv$ is the SAME accepted last-cell quadratic, restored from
saved log/root enclosures and possible-active lists. No logarithm, Gamma
or prime comparison is regenerated. Exact affine quadratic minimization
checks the saved whole-cell minimum and all positive subcell margins;
the outward support strips remain covered. The full reflected partition
has25 unchanged plus128 new cells.

With the original unphased coefficients, let
$g(s)=2|\sum_{2049}^{32768}u_n\cos(\pi ns)|^2$ and let $F$ be its
nonconstant sine primitive. Its two old endpoints are recovered by

$$
F(l)=a_0(1-l)-E_{\rm last},\qquad F(1)=0.
$$

The harmonic intermediate was not archived. The program restores it once
with the existing FLINT convolution/autocorrelation helper, retaining both
sum and difference frequencies, and evaluates only127 missing interior
primitives. It reuses $a_0$ rather than taking a new norm contraction.
All interval coefficients and the accepted anchored trigonometric helper
are retained at1024 working/640 outward bits. Negative original mass
intervals stop before serialization; a valid zero-crossing interval is
intersected with mathematical nonnegativity. The128 new masses must be
consistent with the saved last-cell total. No source/action, inverse,
solve, Gram or far supplier, or old cell integral is replayed.

The same153 masses and weights give the following approximate values:

| Quantity | Value |
|---|---:|
| $a_0=\|u_J\|^2$ | $1344.75921630745$ |
| $E_{\rm last}$ | $1134.58229325125$ |
| $a_1=\int\widetilde W|u_J|^2$ | $10681.5778952432$ |
| $a_2=\int\widetilde W^2|u_J|^2$ | $1164563.05493539$ |
| $m=\|\widetilde W\|_\infty$ | $554.477544320329$ |
| $mU_F$ | $25954.2272562998$ |
| Centered-Cauchy complete upper | $50854.0911649295$ |
| $47h_{\rm low}$ | $50553.3119129287$ |

The centered upper includes
$2\sqrt{U_F}\sqrt{a_2-a_1^2/a_0}$, using unweighted $u_F\perp u_J$.
It remains INCONCLUSIVE, with no actual ineligibility or sign inferred.
The fixed refinement takes5.7965seconds; its128 new mass intervals and
moment bounds are in
[full-gamma-prime-potential-refined.json](full-gamma-prime-potential-refined.json).
There is no refinement/cell/precision/shift scan or retry.

[full_gamma_prime_potential_resolvent.py](full_gamma_prime_potential_resolvent.py)
uses those SAME153 saved masses in a different standard weighted-Cauchy
upper. Its auxiliary price $\lambda=640$ is fixed from the inherited
uniform bound632 before evaluating any resolvent moment. The actual
shift stays $1/48$. Since $\lambda>m$, define

$$
R_\lambda=\int_I\frac{|u_J|^2}{\lambda-\widetilde W}
 =\sum_C\frac{E_C}{\lambda-\widetilde W_C}.
$$

Orthogonality gives $\langle u_J,u_J+u_F\rangle=a_0$. Weighted Cauchy
with the positive multiplier $\lambda-\widetilde W$ therefore yields

$$
a_0^2\le R_\lambda
 \left(\lambda\|u_J+u_F\|^2-\int\widetilde W|u_J+u_F|^2\right).
$$

Using the complete far cap, the single joint upper is

$$
p_v\le\int\widetilde W|u|^2
 \le\lambda(a_0+U_F)-a_0^2/R_\lambda.
$$

The directed implementation uses
$\lambda(a_0^{\rm up}+U_F)-(a_0^{\rm low})^2/R_\lambda^{\rm up}$.
It retains the far self-energy and cross contribution in one common
quadratic form. It assumes neither weighted orthogonality nor jointly
attained separate extrema. Only one scalar sum is new; no coefficient,
primitive or energy is recomputed. Approximately,

$$
R_{640}=2.14857834484827,\qquad U_{\rm res}=48941.0175472545.
$$

Exact exported endpoints certify

$$
\boxed{p_v\le U_{\rm res}<47h_{\rm low}\le47h.}
$$

Thus this same-source fixed1/48 directional budget is conditionally
admitted. The scalar application takes0.0144seconds; its bounds, fixed
price and input/source invocation provenance are in
[full-gamma-prime-potential-resolvent.json](full-gamma-prime-potential-resolvent.json).
No auxiliary price or parameter scan is used.

Both applications require the inherited actual common closed domains,
legal finite-$E$ shifts, original unitary mass/phase identification,
SAME short $\Pi$ and sound full-source/far archives. The earlier old-W
obstruction is unchanged. Budget admission supplies no actual sufficient
score, good-block contribution, retained/cofinal positivity, Robin tail
or RH. These are parameter-specific applications of existing quadratic,
finite-Fourier and weighted-Cauchy facts, with no generic novelty or Lean
verification claim.

Each command reads the published pinned inputs and requires a fresh
output path. The second command consumes the saved153 masses directly.

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_prime_potential_refined.py --out build/c27-prime-potential-refined.json
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_prime_potential_resolvent.py --out build/c27-prime-potential-resolvent.json
```

## Transport the comparison before solving another trial

[full_gamma_shifted_projection_upper.py](full_gamma_shifted_projection_upper.py)
consumes the accepted two-source projection upper $U_2$, original
$r_V=w_0^*R_{\rm res}V$, true positive $h=w_0^*\mathcal Mw_0$
enclosure and the conditional same-source1/48 resolvent budget.
The existing rank-one inverse and inverse-projection facts give

$$
V^*(Q-R_{\rm res}^*\mathcal M_\tau^{-1}R_{\rm res})V
\preceq U_{2,\tau}=U_2+\tau r_V^*r_V/h.
$$

Here $\mathcal M_\tau=\mathcal M+\tau uu^*/[h(1-\tau)]$ and
$u=\mathcal Mw_0$ is the SAME full unscaled source. Every admitted
inverse-majorant/residual-trial lower construction for this comparison
is dominated by the left side. The inequality requires the inherited
actual common inverse forms/domains, legal $E$ elimination, unitary
identification, SAME short and sound common source/lift/full-far
archives. A finite source projection does not truncate the exterior.

In particular, transporting
$\overline W_\tau=\overline W-\tau w_0w_0^*/h$ along with
$w_\tau=(1-\tau)w_0$ preserves $\mathcal B$ and $\mathcal K$,
while its base gains $\tau r^*r/h$. This lower dominates the
held-majorant shifted-row construction. Both remain below the exact
shifted comparison and hence below the displayed compressed ceiling.

The new packet
[full-gamma-shifted-projection-upper.json.zlib](full-gamma-shifted-projection-upper.json.zlib)
uses1024-bit arithmetic and the existing640-bit entry grid, with
11584-bit determinant export. Its one new determinant satisfies

$$
\det U_{2,1/48}<-4.93328\times10^{-1362}<0.
$$

Under the stated premises this excludes whole-retained positivity
certification by the specified fixed $\mathcal M_{1/48}$ comparison,
including both transported and held-majorant lower constructions.
It is not actual Weil negativity. The run takes0.0286seconds, with
no source/action, energy, projection, inverse, solve, Gram or far
production. A new seventeen-dimensional trial is not needed and
was not executed. No parameter/source/precision scan or retry occurs.
No actual retained/cofinal positivity, Robin, RH, Lean or generic
novelty conclusion follows.

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_shifted_projection_upper.py --out build/c27-shifted-projection-upper.json.zlib
```

### The entire supplied eligibility interval

[full_gamma_shifted_projection_interval.py](full_gamma_shifted_projection_interval.py)
uses only exact rational scalar arithmetic on the two archived
zero/1/48 determinant enclosures and the saved resolvent upper
$U_{\rm res}\ge p_v$. With the original directed $h_{\rm low}$,

$$
\tau_{\rm cap}=\frac{h_{\rm low}}{h_{\rm low}+U_{\rm res}}
\approx0.0215048994247574
$$

supplies eligibility for every $0<\tau\le\tau_{\rm cap}$ under
the SAME actual common-form/short/source premises. The resolvent
bound concerns $p_v$ independently of the shift; its original target
at1/48 does not restrict this implication to that single point.

A rank-one determinant is affine even for a singular base. Thus
$d(\tau)=\det(U_2+\tau r_V^*r_V/h)$ obeys

$$
d(\tau_{\rm cap})=(1-48\tau_{\rm cap})d(0)
 +48\tau_{\rm cap}d(1/48).
$$

Because $\tau_{\rm cap}>1/48$, its directed upper uses the lower
endpoint of $d(0)$ and the upper endpoint of $d(1/48)$. The maximum
of the two endpoint uppers then bounds the full affine segment.
The new packet
[full-gamma-shifted-projection-interval.json](full-gamma-shifted-projection-interval.json)
contains the exact cap and determinant bounds and certifies

$$
\sup_{0\le\tau\le\tau_{\rm cap}}d(\tau)
<-4.36494\times10^{-1362}<0.
$$

The scalar run takes0.0145seconds, without any new determinant,
inverse, solve, projection, source/action, energy, Gram or far supplier.
Consequently the specified fixed-source $\mathcal M_\tau$ comparison
cannot certify whole-retained positivity anywhere in this supplied
interval. Numerical refinement of a residual trial inside the same
comparison cannot alter that ceiling.

This is a sufficient budget interval, not the largest actual capacity.
It does not exclude larger actually eligible shifts, other update
directions or changed comparison/support. All inherited actual
inverse-form/domain, legal $E$, unitary/SAME-short/source/lift/full-far
and sound-archive premises remain conditions. Actual retained/cofinal
positivity, Robin and RH remain unresolved; no Lean or generic novelty
is claimed.

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_shifted_projection_interval.py --out build/c27-shifted-projection-interval.json
```

## Full-band orthogonality for the same-source surplus budget

[full_gamma_band_surplus_budget.py](full_gamma_band_surplus_budget.py)
uses the SAME original384-bit source,153-cell reflected refined weight
$W\ge1/v$, complete far cap and saved moments. In the original
unphased orthonormal cosine basis, $u_J$ occupies2049..32768 and the
entire $u_F$ is perpendicular to every mode0..32768. With $P_N$ the
projection onto that complete band, the standard Bessel/Cauchy bound is

$$
p_v\le a_1+2\sqrt{U_F}\sqrt{a_2-\|P_N(Wu_J)\|^2}+mU_F.
$$

For the reflected coordinate $s$, define
$c_l=\int_0^1W(s)\cos(\pi ls)\,ds$. Rational step jumps give every
needed $c_l$, and the normalized output coefficients are

$$
\beta_0=\sqrt2\sum_{n\in J}u_nc_n,\qquad
\beta_k=\sum_{n\in J}u_n(c_{|n-k|}+c_{n+k}),\quad1\le k\le32768.
$$

One native ball polynomial product implements the symmetric Laurent
convolution, keeping both sum and difference frequencies. The input
helper requires positive source modes; it does not claim constant-source
support. Its zero OUTPUT mode has the separate $1/\sqrt2$ factor.
The exact finite multiplication identity uses only weight coefficients
through65536; it does not approximate the step weight by a finite series.
This reuses Olver–Townsend,
[arXiv:1202.1347v2](https://arxiv.org/pdf/1202.1347v2), §2.2, (2.7),
and Johansson, [arXiv:1611.02831v1](https://arxiv.org/pdf/1611.02831v1),
§6.1. No continuous-coefficient uniform approximation or unverified
floating-point FFT error bound is imported.

The prescribed sole1024/640-bit run takes9.168seconds and supplies:

| Quantity | Approximate enclosure value |
|---|---:|
| $\|P_N(Wu_J)\|^2$ |1147558.75453915025|
| $\|(I-P_N)Wu_J\|^2$ upper |17004.3003962366169|
| Complete cross upper |1784.31468694854120|
| Complete $p_v$ upper |38420.1198384914721|
| Supplied eligibility cap |0.0272333893565792018|

The Bessel subtraction uses the projection LOWER and saved $a_2$
UPPER. The final bound includes the entire $mU_F$ and cross contribution.
No old source-cell mass, coefficient source/action, operator, potential,
inverse, solve, Gram or far producer is replayed. All65537 new weight
coefficients and32769 projected coefficients are retained as finite
ordered dyadic intervals, enabling subsequent consumers to reuse them.

The same run reuses the old zero and1/48 determinant enclosures to
enclose their affine zero near0.026662647029077637 and the necessary
capacity target near39265.56779425506. Its primary success criterion
is the STRICT complete-budget inequality against that target's LOWER
endpoint. The new upper passes. Affine propagation at the sufficient
new cap gives a positive determinant near4.83012948603e-1363, with
no new determinant/matrix calculation. This escapes the preceding
ceiling exclusion only; a positive upper ceiling is not a positive
lower or an actual retained-positivity certificate. No maximal actual
capacity is claimed. The previous exclusion on its smaller supplied
interval remains valid.

Accepted input byte pins and correspondence are exported in
[full-gamma-band-surplus-budget.json.zlib](full-gamma-band-surplus-budget.json.zlib),
SHA256 `04fb6eb533b8448ade52ad4032faaaf079ad6603f9998d272a36828c57be44b1`.
Its source SHA256 is
`4ce09199999cf5730623507c38be4e5301472ee0d4b882947f3fab76c8e8245a`.
Source/helper hashes are invocation guards, not replacement data-cache
keys. All coefficient exports reject nonfinite or unordered intervals;
an inconsistent negative Bessel residual stops rather than being hidden.
Failure of the strict target is inconclusive and ends this fixed recipe,
without a cell, source, parameter or precision scan.

The actual common inverse forms/domains, legal $E$ elimination, same
closed short, unitary source/phase normalization, full-band orthogonality
and sound entire-far/cache identifications remain assumptions. These are
applications of established methods. Actual sufficient score/good-block,
retained/cofinal positivity, Robin, RH and Lean verification are unresolved.

To reproduce into a fresh output, rather than overwrite the archive:

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_band_surplus_budget.py --out build/c27-band-surplus-budget.json.zlib
```

## Reuse the fixed trial at the supplied band endpoint

[full_gamma_scalar_transport.py](full_gamma_scalar_transport.py) consumes
the exact supplied band cap and BOTH true-$h$ endpoints. It reuses the
same exact17trial $x$, old signed-trial LOWER, gradient energy, first
gradient coordinate, mixed source scalar and curvature. The old
held-majorant program supplies these scalar enclosures only; its
verdict and far-penalty constants are not transported.

Write $\Sigma_0=\Sigma_*+\delta H_{\rm lift}$, $a_0<0$,
$g_0=((b_0-cx)/a_0,x)$, $m_0=rg_0$, $z_0=g_0^*Gg_0$ and
$\mu=\tau/h$. The standard rank-one identities for this SAME real
trial give

$$
\eta=1+\mu r_0^2/a_0,\quad a_\tau=a_0+\mu r_0^2,\quad
\theta=\mu r_0m_0/a_\tau,\quad g_\tau=g_0-\theta e_0,
$$

$$
q_\tau(x)=q_0(x)-\mu m_0^2/\eta,\qquad
z_\tau=z_0-2\theta(e_0^*Gg_0)+\theta^2G_{00}.
$$

Only the missing18term cross $e_0^*Gg_0$ is newly contracted. The raw
complete paid $G$ is used; the energy-scaled joint Gram is not substituted.
For $\kappa=33/32$, the NEW coefficient and directed exclusion lower are

$$
A_\tau=\kappa-\kappa^2G_{00}/(-a_\tau)>0,\qquad
\ell_\tau=\mathrm{charge}_{\rm lo}+q_{0,\rm lo}
 -(\mu m_0^2/\eta)_{\rm up}-(\kappa^2z_\tau/A_\tau)_{\rm up}.
$$

The sole fixed1024/640-bit application takes0.09882seconds. Its
[saved packet](full-gamma-scalar-transport.json), SHA256
`c87498beff259c4a22ba319aa32af4975649c35a50990ba27bb2324a43ffa757`,
records the nine accepted input pins, all14 invoked source/helper hashes,
true-$h$ bounds, exact cap and new scalar enclosures:

| Quantity | Approximate enclosure value |
|---|---:|
| $\tau$ |0.0272333893565792018|
| $\eta$ |0.379630426033288183|
| $a_\tau$ |$-4.47097042972572520\times10^{-32}$|
| $A_\tau$ |0.362304559585549403|
| $z_\tau$ |55.6410065569858228|
| Signed-trial decrement upper |60.8094277285621680|
| New completed-square penalty upper |163.323659119699822|
| Exclusion lower $\ell_\tau$ |$-177.935387035771641$|

The pivot and far-coefficient guards pass, but the exclusion lower is
negative. The result is INCONCLUSIVE: it proves neither feasibility nor
actual positivity or negativity. If a strict positive endpoint lower
were obtained, the unchanged $B,K$ and fixed-$D$ Loewner monotonicity
would exclude the SAME transported family throughout the supplied
interval; that implication is not triggered here.

No new trial, matrix, determinant, inverse, solve, Fourier/potential,
source/action, Gram or old energy/far producer is evaluated. No parameter
or precision scan follows this inconclusive result. A sharper exclusion
would need a justified further estimate; a sufficient positivity
certificate still needs the SAME actual Schur complement and good-block
inverse cost, which the projected18row does not supply. The actual
common forms/domains, legal $E$, closed short, unitary/source/lift and
sound full-far/cache premises remain conditions. These are established
rank-one and square-completion methods applied to saved data, not a new
generic theorem, maximal capacity, Lean, Robin-tail or RH result.

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_scalar_transport.py --out build/c27-scalar-transport.json
```


## Assemble the complete retained finite row from cached arrays

[full_gamma_finite_retained_row.py](full_gamma_finite_retained_row.py)
assembles the original real, unphased384-bit row0 finite contribution
on retained modes $0,\ldots,1280$:

$$
b_J=c^*-(A_\rho^{-1}F_\phi^*)^*\Beta.
$$

It reads `common.column` (1281-by1),
`actual.channel_inverse_actual_column` (19-by1), and
`boundary.finite_beta` (19-by1281). One thin channel product supplies
all1281 coefficients. One product with the saved1263-by18 inverse
action checks the SAME exact old lift $V=[-X;I]$ in its saved order:

$$
b_JV=(b_J)_\ell-(b_J)_gX.
$$

The sole1024-working/640-export-bit assembly takes0.22688seconds.
All1281 finite ordered entries are exported, and all18 projected
coordinates overlap the archived actual finite lift row. Source identity
and distributivity, rather than interval overlap alone, identify the
common mathematical row. Packet SHA256:
`479b8da2bb1d490d0f4adf285a67797d733ec05d4dab33d88b6175946a89f28a`.

The [saved finite-row packet](full-gamma-finite-retained-row.json.zlib)
contains the1281 coefficients,18 lift-check coefficients, six fixed input
pins, invoked source/helper hashes, complete coordinate order and
inherited conditions. Fresh output, finite ordered exports and unchanged
inputs/helpers are required. No old producer, source/adjoint/channel
projection, inverse, solve, Gram, far or Fourier calculation is replayed;
no parameter, source or precision scan runs.

This is complete only in retained coordinates. The actual full
$b=\mathcal B$ includes an additional SAME-source far mixed row. The
full unweighted $r=w^*R$ and a certified1263good-block inverse metric
remain missing. The cached17inverse belongs inside the18-core; neither
it nor1263-by18 saved actions gives an inverse for a new good-block row.
The [owning application note](https://github.com/the-omega-institute/trureturing/blob/dev/Library/Weil/connesconsanimoscovici2026spectral.md)
states the standard sufficient two-row bordered-inertia interface and
its remaining complete-row/metric inputs. Those methods are reused,
not claimed as new results. The assembly supplies no sufficient score,
actual retained/cofinal positivity, Robin-tail bound, RH or Lean result.

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_finite_retained_row.py --out build/c27-finite-retained-row.json.zlib
```

## Stop the fixed-charge true-old-base recipe before filling its metric

[full_gamma_old_base_stop.py](full_gamma_old_base_stop.py) uses the TRUE
old18-core $\Sigma_0$, the original unweighted source row $\rho=rV$,
the finite row $\beta_J=b_JV$ and the cached exact17-coordinate trial
$x$. It does not substitute the earlier compression upper
$\Sigma_{\max}$. Its cached trial lower is valid for this exact $x$:
the producer's positive variational lower was strictly positive, so
the `max(0,trial)` clipping was inactive.

Partition $\Sigma_0=\left(\begin{smallmatrix}a&c\\c^*&F\end{smallmatrix}\right)$,
and reuse its certified inertia$(17+,1-)$. Choose exact rational values

$$
t_+=\tau_{\rm cap}/h_{\rm lower},\qquad
\overline D=\operatorname{up}(\texttt{refined_common_charge_upper}).
$$

The first coefficient is an UPPER COMPARISON, not an eligibility claim.
Every actual $t=\tau/h$ in the supplied band satisfies $0\le t\le t_+$.
The only missing contractions restore

$$
g_0=\binom{(\beta_{J0}-cx)/a}{x},\qquad m=\rho g_0.
$$

With $\eta=1+t_+\rho_0^2/a>0$, the established fixed-trial identity gives

$$
q_{\mathrm{trial},t_+}
=q_{\mathrm{trial},0}-t_+m^2/\eta,
\qquad
q_{\mathrm{trial},t_+}\ge
q_{0,\mathrm{lo}}-\operatorname{up}(t_+m^2/\eta).
$$

The sole1024-working/640-export-bit run takes0.10488seconds:

| Quantity | Approximate enclosure value |
|---|---:|
| $t_+$ |$2.5319197720691323970\times10^{-5}$|
| $a+t_+\rho_0^2$ |$-4.47097042972572520\times10^{-32}$|
| $\eta$ |0.379630426033288183|
| $q_{0,\mathrm{lo}}$ |$-120.742434607587509$|
| Transport loss upper |60.8094277285621680|
| Transported trial lower |$-181.551862336149677$|
| Fixed $\overline D$ |183.800677119682971|
| Stopping margin lower |2.24881478353329417|

The negative updated first pivot preserves exactly one negative core
direction. The strict positive stopping margin makes the center-only
two-row border positive definite; adding ANY PSD joint good-block
metric preserves it. The existing inertia criterion therefore excludes
PSD of the auxiliary target with zero PROJECTED far response. Target
order extends the exclusion to all $0\le t\le t_+$ and
$D\ge\overline D$ for that same recipe.

The [saved directed packet](full-gamma-old-base-stop.json.zlib), SHA256
`685d3535cfbc21f6ad3d03b04126d5250ff76570265fbf7f315ca741a6b45af3`,
pins seven inputs, nine invoked source/helper files and seven cached
producer sources. Only17+18 contraction terms are new. No producer,
inverse, solve, trial quadratic, profile, Gram, source or far action is
replayed, and no parameter/source/precision scan follows.

An enlarged energy/Gram UPPER-only far-response relaxation that permits
zero projected response cannot certify this fixed-charge recipe
uniformly. Tighter upper envelopes that still permit zero do not remove
this obstruction. Zero response is not identified with the actual
fixed source. A source-specific signed/correlated far-charge estimate,
a smaller valid denominator, or a stronger TRUE baseline can escape it.
The actual source far completion, $D<\overline D$, retained/cofinal
positivity, Robin tail, RH and Lean remain unresolved. This is a new
cached-data application of established variational/inertia identities,
not a new general method or a rerun of the unshifted exclusion.

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_old_base_stop.py --out build/c27-old-base-stop.json.zlib
```


## Signed same-source infinite Hankel tail energy

[full_gamma_hankel_tail_energy.py](full_gamma_hankel_tail_energy.py)
uses the accepted exact384-bit row0 w and enclosed q to bound the SAME
actual action $u_F=G_{FJ}w-G_{FE}q+H_{FE}(q-p)$. It combines the signed
principal sources before taking a norm. The standard harmonic
partial-fraction tail Gram is evaluated by one FLINT linear Toeplitz
convolution; only its diagonal is padded by the integral bound.
The complete derivation and parameter correspondence are in the
[owning application exposition](https://github.com/the-omega-institute/trureturing/blob/dev/Library/Weil/connesconsanimoscovici2026spectral.md#retain-signed-source-correlations-in-the-infinite-principal-tail).

The [directed result](full-gamma-hankel-tail-energy.json.zlib), SHA256
`d09a5ab254d164cfd40982b31e1a155e31d6288b705d916bf248236a0d527606`,
records one1024-working/640-export-bit application taking0.52936seconds:

| Quantity | Previous upper | New upper |
|---|---:|---:|
| Complete weighted far energy |14.1648613043077766|9.16483973939069465|
| Positive source-valid fee |183.800677119682971|178.644404880862230|
| Supplied sufficient parameter cap |0.0272333893565792018|0.0358711862066475850|

The norm factor is1/2, so the principal energy factor is1/4. Saved
smooth remainder bounds already include gamma and are UNWEIGHTED;
they are multiplied once by sqrt(Dfar). The saved mixed-H norm is
already WEIGHTED and remains paid in full. Rebuild the charge constant
from the original middle upper, h lower and channel-energy lower;
subtracting an old far-energy upper from an old charge upper is invalid.
A finite positive fee is checked before use. The saved full-band Bessel
projection, finite weighted moment and maximum reflected weight are
reused unchanged for the cap.

The packet pins nine accepted inputs, nine invoked source/helper files
and nine cached producers at entry and exit. It saves the preselected
convolution coefficients0,767,31487 for a thin result check. No source,
inverse/solve, far-mode or old producer is regenerated; no new Fourier
potential projection, determinant or parameter/source/precision scan
is performed. The previous finite w*Jw convolution is a different
quadratic and is not replayed. The prescribed recipe stops on failed
or nonfinite guards or no strict energy improvement.

These are conditional actual-source energy/fee/eligibility estimates.
They preserve all inherited actual-form/domain, inverse/legal-E,
SAME-short, unitary/source/lift/comparison and sound full-far/cache
premises. The smaller fee is outside the previous fixed-fee exclusion,
but gives no actual FULL border sign. Complete signed far response,
joint good-block inverse metric, actual retained/cofinal positivity,
Robin, RH and Lean verification remain unresolved. Classical tail
identities and existing FLINT/Bessel methods are reused, with no
new generic theorem or algorithm claim.

For an explicitly requested reproduction, choose a fresh output path:

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_hankel_tail_energy.py --out /tmp/c27-hankel-tail-energy-reproduction.json.zlib
```

## Fixed negative-channel envelope capacity

[full_gamma_negative_channel_capacity.py](full_gamma_negative_channel_capacity.py)
tests one fixed full-center UPPER-envelope package before supplying
an actual signed far moment. It reuses the TRUE old18-core, its exact17
reference inverse and certified Neumann bound, the complete paid far
Gram and saved trace, the signed-tail energy/fee/cap, and original
unphased384-bit row0. The legal coefficient is cap/h_upper; the old
energy-cap lower used to normalize the trace is not an actual
source-energy lower. The
[owning application exposition](https://github.com/the-omega-institute/trureturing/blob/dev/Library/Weil/connesconsanimoscovici2026spectral.md#a-fixed-negative-channel-envelope-before-a-signed-far-functional)
gives the block inverse decomposition and the directed endpoint rules.

Three17-dimensional reference forms select one tail-independent64-bit
dyadic alpha and price its positive channel. For that fixed direction,
write the actual finite negative-channel response as lambda and the
unknown far response as xi. The resulting upper envelope is
Tbar-(lambda+xi)^2/delta, with |xi|<=Cminus. Its directed stopping
floor uses delta_lower; the future sufficient response threshold uses
delta_upper.

The [saved packet](full-gamma-negative-channel-capacity.json.zlib),
SHA256 `1487736c15e524c3b86859c3e427c85ccf3512659990597bb8199e882688b21e`,
records one1024-working/640-export-bit run taking0.04108seconds:

| Quantity | Approximate directed value |
|---|---:|
| Fixed alpha |0.0379872985139444681|
| Legal t0 |0.0000333498575645498988|
| Frozen positive envelope Tbar |911.633113280467383|
| Finite negative-channel response lambda |9.64684697797322708e-15|
| Allowed far-response radius Cminus |5.23550628472335536e-16|
| Sufficient response magnitude upper |1.03616862656981840e-14|
| Package stopping floor lower |33.3488964810889399|

The strictly positive floor excludes strict negativity of this ONE
fixed upper-envelope package for every far response allowed by its
energy/Gram bounds. An actual signed xi alone therefore cannot make
this package certify negativity; its prospective signed-tail supplier
is not executed. Tighter coherent positive-channel pricing, another
direction or another legal recipe is outside this exclusion. Even the
most favorable separate endpoints were used for the stop, with no
assertion that they share an actual realization.

The packet pins ten accepted inputs, invoked local helpers and ten
cached producer sources at entry and exit. No old producer or passed
suite, inverse/solve, trace/Gram/source/far supplier or determinant is
replayed, and there is no parameter/source/precision scan. Source and
result reviews use new tiny implementation examples and thin scalar
consumption only. This is a conditional cached-data application of
standard block inverse, Loewner and norm inequalities, with no generic
novelty or Lean claim. Actual full-border sign, the additional good-block
metric, retained/cofinal positivity, Robin tail and RH remain unresolved.

For an explicitly requested reproduction, use a fresh output path:

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_negative_channel_capacity.py --out /tmp/c27-negative-channel-capacity-reproduction.json.zlib
```

## One coherent same-source tail direction

[full_gamma_coherent_mixed_capacity.py](full_gamma_coherent_mixed_capacity.py)
retains the original unphased384-bit row0, fixed alpha, legal t0,
positive fee and source-energy cap. On the SAME true old18-core,
let Q_R=T R T*−e0 e0*/delta, y=alpha rho+betaJ and vR=Q_R y*.
The actual signed functional theta=f vR combines the positive and
negative cross terms before either is bounded separately. Existing
relative Loewner and PSD bounds yield the center upper Clin+2theta,
before the additional PSD good-block metric. R remains a reference,
not the actual inverse. The negative tail self-square is discarded.

The [coherent packet](full-gamma-coherent-mixed-capacity.json.zlib),
SHA256 `b082eaad102c5f7ca38d1c49c19970952288ce630198a4ea76b48b96d8d23b65`,
freezes the packet's exact rational Clin (approximately25.528336100647874028)
and the sufficient signed condition theta<−Clin/2 (threshold approximately−12.764168050323937014).
It reads the cached RAW reference energy,
not the already-inflated positive-channel energy. One new17-reference
vector action and one18-dimensional full paid Gram form take0.04560
seconds at1024working/640export bits. All old producers, identities,
tests, inverse/solve, trace and source suppliers are reused.

The direct Gram interval spans approximately ±4.4981e24 and gives only
the response-radius upper6.62128273610463133e12. Its negative endpoint
is not a negative PSD norm; the inconclusive package result does not
establish a large actual norm, improvement, sign or favorable response.
Saved coefficient intervals enclose ONE exact analytic direction;
they are not independently selectable coefficients or free midpoints.

## Cached PSD trace bound for the coherent direction

[full_gamma_coherent_split_capacity.py](full_gamma_coherent_split_capacity.py)
consumes only the two accepted packets' scalars and metadata. With
the SAME complete paid PSD majorant G, the existing normalized trace
upper Tfar bounds tr(Qplus G). The PSD operator<=trace inequality and
R<=(1+eta)S^-1 give the uniform direction-norm upper

```text
Nsplit = UP(sqrt(r_UP*(1+eta_UP)*Tfar_UP)
            + abs(lambda)_UP/delta_LOWER*sqrt(G00_UP)).
```

Tfar is already normalized; there is no additional energy division.
The norm uses abs(lambda) UPPER and delta LOWER, opposite to the
negative finite credit inside Clin. Take the minimum with the saved
direct norm upper on this SAME analytic direction, then freeze
Ctheta=UP640(kappa*sqrt(Ecap)*Ncap). This is a uniform operator bound;
a particular-source |theta| bound cannot be divided to recover a norm.
The complete derivation and parameter map are in the
[owning application exposition](https://github.com/the-omega-institute/trureturing/blob/dev/Library/Weil/connesconsanimoscovici2026spectral.md#bound-the-same-direction-using-the-cached-positive-trace).

The [scalar-only packet](full-gamma-coherent-split-capacity.json.zlib)
has SHA256 `6e4232d6079a4279511a67c13c6501d51a3b4fdaf410ccfee957adce29e1d7bc`
and records one fixed1024-working/640-export-bit evaluation taking0.01185seconds:

| Quantity | Approximate directed value |
|---|---:|
| Positive-direction norm upper |15.3617882831078431|
| Negative-direction norm upper |13.7365253063170975|
| Selected direction norm upper |29.0983135894249406|
| Response-radius upper |90.8435754570087271|
| Unchanged frozen Clin |25.5283361006478740|
| Package floor lower |−156.158814813369580|
| Uniform center upper |207.215487014665328|
| Unchanged sufficient signed threshold |−12.7641680503239370|

The radius is smaller than the direct interval's bound, but the negative
floor means only that this fixed linear package is not excluded. No
actual signed theta or favorable realization is supplied. All source,
lift, full-far allowance, trace and comparison assumptions remain
conditions. A future same-source signed interval must have an upper
below the frozen threshold and retain the analytic direction relations;
the additional good-block metric would still have to be paid.

The scalar adapter guards both packets, original input/helper/producer
bytes and current invoked sources before and after use. No matrix or
vector is extracted; no old suite/producer, new reference/Gram/trace,
inverse/solve or source/far supplier runs. There is no source, parameter
or precision scan. Nonexclusion alone does not justify an expensive
supplier; its factored or short-lift interface and complete error
contract must first be established. Existing classical block inverse,
Loewner, PSD trace and Cauchy bounds are applied without a new generic
theorem, algorithm or Lean claim. Actual full-border sign, good-block
metric, retained/cofinal positivity, Robin tail and RH remain unresolved.

For an explicitly requested reproduction, use a fresh output path:

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_coherent_split_capacity.py --out /tmp/c27-coherent-split-capacity-reproduction.json.zlib
```

## Same-source signed-jet interface

[full_gamma_signed_jet_interface.py](full_gamma_signed_jet_interface.py)
consumes the saved full56 features F_lift=P0 Dphys V and the ONE
coherent direction vR. The physical residual Dphys=[I;-Y] includes
modes0..2048; the far synthesis starts at32769 with scale M=32768.
One cached112x18 action encloses h=F_lift vR, in order
A0..A55,B0..B55. All source/lift coefficient intervals propagate;
there is no independently selected midpoint direction.

The original geometric and exact-midpoint rounding components imply
the forward bound ||BvR-j_h||_W<=sqrt(2eps_complete)*||VvR||. Here r2
is the original UNWEIGHTED squared Frobenius error. eps_complete
already pays the original far weight and physical residual amplification;
the original allowance metadata grid is512, while matrix grids are640.
A Gram-error estimate alone would not supply this forward bound. The
genuine full moment Gram controls the feature error; the masked core
is not assumed PSD.

Eighteen diagonal lift-Gram terms bound ||VvR|| by
LV=UP sum_j abs(vR_j)_UP sqrt(liftGram_jj_UP). Freeze
eps_b=UP640(sqrt(2eps_complete)*LV). For the SAME principal proxy u0,
reuse the independent principal norm upper P and smooth+mixed-H norm
error eps_u. With theta=kappa<uF,BvR> and mu=kappa<u0,j_h>, telescoping
gives error E_mu=UP640[kappa*(eps_u*Ncap+P*eps_b)]. The proxy uses
minus q on1281..2048 and the original w on2049..32768; a q remainder
is not renamed a p remainder.

The [saved interface](full-gamma-signed-jet-interface.json.zlib), SHA256
`154949d8ab3e1f2875de531aa63a9c5a3c3d7d133640a5ce8b1d23c98d67565a`,
records one fixed1024-working/640-export-bit evaluation taking0.28065seconds:

| Quantity | Approximate directed value |
|---|---:|
| Lift triangle bound LV |4.48767188714705487e67|
| Inherited complete allowance |2.07732485158467804e-135|
| Forward jet error eps_b |2.89259947122008663|
| Principal-proxy norm error eps_u |0.478756042962487477|
| Signed-functional error E_mu |21.9687713914200247|
| Unchanged Clin |25.5283361006478740|
| Sufficient proxy threshold |−34.7329394417439617|

The112 coefficient intervals have maximum absolute endpoint
approximately1038.04984493700085 and maximum absolute width
approximately5.11418e-25. The exact future sufficient condition is
UP(Clin+2*(I_mu+E_mu))<0, or in exact arithmetic
sup(I_mu)<−Clin/2−E_mu. The exact rational threshold is saved; decimals
are approximate. No actual mu, theta, center or full-border sign is
supplied. The extra PSD good-block metric remains unpaid.

This adapter reuses the existing component, norm and source certificates
and passed suites. It guards source/helper/producer/input bytes at entry
and exit and writes one exclusive fresh archive. No reference, Gram,
trace, inverse/solve, source/far supplier, old producer or scan runs.
It applies standard triangle and Cauchy inequalities without a new
generic theorem, algorithm or Lean claim. The common actual-form/domain,
kernel, unphased source384/row0, short/lift, original weight and complete
component-bound premises remain conditions. A priced threshold alone
does not justify an expensive supplier. Actual signed response,
good-block metric, full-border sign, retained/cofinal positivity,
Robin tail and RH remain unresolved. The full derivation is in the
[owning application exposition](https://github.com/the-omega-institute/trureturing/blob/dev/Library/Weil/connesconsanimoscovici2026spectral.md#use-the-full-short-residual-jet-for-one-signed-functional).

For an explicitly requested reproduction, use a fresh output path:

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_signed_jet_interface.py --out /tmp/c27-signed-jet-interface-reproduction.json.zlib
```

## Same-source signed proxy with a bounded oscillatory tail

[full_gamma_signed_mu.py](full_gamma_signed_mu.py) applies the existing
FLINT reciprocal-sum convolution, iterative polynomial evaluation and
real Hurwitz-zeta supplier to the original384-bit source and all112
saved full56 jet coefficients. The source uses minus q on1281..2048
with q=`D_inverse_difference`, and the original w on2049..32768. The
normalization is M32768, the finite output block is32769..65536, and
the analytic tail starts at65537. No source or direction is reselected.

The scalar identity, source-power expansion, and complete source/symbol
error telescope are in the owning
[application exposition](https://github.com/the-omega-institute/trureturing/blob/dev/Library/Weil/connesconsanimoscovici2026spectral.md#enclose-the-same-source-signed-proxy-with-an-abel-tail-allowance).
For F0=a_r[xB(x²)-C(x)A(x²)] and F1=a_r A(x²), only the F0 real Hurwitz
tail is included in the computed center. Classical Abel summation
bounds each omitted unit-phase tail by
2*M^(q-1)/(|1-xi|*T^q). Positive phase gaps are certified, and the
resulting Rosc is added to the original R_mu. The bound uses the SAME
approximate sourceF1; no actual-source bound is substituted into a
different telescope. All zero-crossing coefficient intervals remain.

A fresh exclusive Cauchy checkpoint saves32768 reciprocal samples
before finite-symbol evaluation, and a second checkpoint saves32768
approximate-jet outputs before infinite-tail evaluation. Finite sine
samples use FLINT real sin_pi(m*2log(u)/log27) directly, retaining
interval uncertainty without rectangular rotation wrapping. The
final archive references both checkpoints by SHA256 and saves the new
real kernel intervals and certified phase-gap lower bounds for reuse.
The finite approximate jet alone is not an actual-operator certificate;
all symbol/forward errors remain paid in the final interval contract.

The saved direct-sine/Abel application uses1024 working bits and640-bit
exports, taking6.96478seconds. Its fixed source-power order is13;
it evaluates186 real Hurwitz kernels,13 positive phase gaps and425984
direct real sine samples. One reciprocal convolution supplies both
the finite functional and the reusable Cauchy checkpoint.

| Quantity | Approximate directed result |
|---|---:|
| Finite approximate functional |−2.266706962449528562|
| Smooth tail center |−2.242939503220959819|
| Algorithmic center |−4.509646465670488381|
| Original source/symbol allowance R_mu |0.001800203553588446|
| Omitted oscillatory allowance Rosc |0.000869186655284127|
| Complete analytic allowance |0.002669390208872574|
| Exported proxy I_mu |[−4.512315855885066397,−4.506977075455910365]|
| Exported actual-theta enclosure |[−26.48108725775422862,17.46179432641325185]|
| Existing Clin+2theta enclosure |[−27.43383853406987275,60.45192487268366729]|

Decimals are approximate; exact grid intervals/rationals are authoritative.
Arb radius construction conservatively inflates the exported intervals.
The proxy interval does not meet the inherited threshold−34.73293944.
Moreover, exact endpoint arithmetic gives Clin+2I_mu>16.5 throughout.
Reducing the nonnegative transfer allowance E_mu alone, even to zero,
cannot make this proxy-only sufficient upper test negative. This does
not determine the actual center's sign. For delta=theta−mu, a sufficient
new signed bound is sup(delta)<−Clin/2−sup(I_mu), approximately
−8.25719097486802665. The inherited absolute delta bound does not
supply that negative upper.

The final [interval packet](full-gamma-signed-mu-direct-sine.json.zlib)
has SHA256 `6475e5f63dba423b2f1ee7e488581f15ef8d53de83a22ea064ac2e393b124d55`.
It references the [Cauchy checkpoint](full-gamma-signed-mu-direct-sine.json.cauchy.zlib),
SHA256 `dbfedf8f55f74b4fa1659ff02e57559fea52e3250379c8f99ed1807ee91245f9`,
and [finite-jet checkpoint](full-gamma-signed-mu-direct-sine.json.finite.zlib),
SHA256 `43a29c6710cae4b4b7ad2fce9433ebba4ca0bdec9dfea968d668b13d3a469bc8`.


All reported intervals retain the original same-source/lift/domain and
component-bound premises. The complete sufficient center test still
includes the inherited E_mu and the frozen Clin; the additional PSD
good-block metric is separate. No positive sufficient upper expression
is interpreted as a positive actual center. No full-border, retained
or cofinal positivity, Robin-tail or RH result, new generic algorithm
or Lean theorem is supplied. Old producers and passed suites are
reused, with no reference/Gram/trace/inverse/solve, source, parameter,
precision or cutoff scan.

For an explicitly requested reproduction, choose a fresh final path
whose corresponding `.finite.zlib` and `.cauchy.zlib` checkpoints
also do not exist:

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_signed_mu.py --out /tmp/c27-signed-mu-abel-reproduction.json.zlib
```

## Same-source signed mixed correction from the saved Gram

[full_gamma_signed_correction.py](full_gamma_signed_correction.py) consumes
the saved proxy interval and original physical Gram kernel, full112 short
jet, inverse profiles and scalar/component bounds. For q=`D_inverse_difference`
and p=`F_inverse_z`, lambda=q−p is supported on1281..2048. The original
source identity is uF=u0+r_sm+Cphys*ell; the target is Bv=Cphys*t, with
cached h=P0*t. The missing signed cross is
nu=kappa*<Cphys*ell,Cphys*t>W, enclosed around kappa*a^T*K0*h,
where a=P0*ell. The original physical `far_matrix_error` pays the cross
error; it is read on its512-bit metadata grid, rather than replaced by
the amplified complete allowance. The saved kernel is not assumed PSD.

Two short calls to the existing signed source-power reducer supply a:
even coefficients of (−1)^n*lambda_n and odd coefficients of
−(−1)^n*lambda_n*s_mid_n, using its negative power ratio. All original
kernel entries and right features are reused. The old primal Hmix bound
and a positive lower bound on its symbol coefficient supply a lambda-l1
upper without another absolute source reduction.

The existing separate geometric/feature forward components also give
a sharper proxy-forward allowance. The true full moment norm is read
on its628-bit grid. The physical symbol-width premise is checked on
all2049 rows; r2 is square-rooted at1024 bits before640-bit rounding,
and its squared interval is retained separately on a1024-bit grid.
The complete allowance for theta around saved I_mu+nu0 is
kappa*(eps_sm*Ncap+P*eps_b_sharp)+kappa*eps_far*Zlambda*sqrt(Kphys)*LV.
This replaces the earlier E_mu; I_mu already includes the proxy's
source, symbol and Abel allowances.

The sole fixed1024/640 application finishes in1.30498seconds. It performs
two short source-power reductions on768 entries and one112-by112 cross.
No old proxy/source, scalar, moment, inverse, reference, trace or Gram
producer is replayed. Its new left features and fine squared remainder
are retained in the [interval packet](full-gamma-signed-correction.json.zlib),
SHA256 `d527c5217314bdd5b21e0f1080ef17b1e29b880954a1bafe77a711a1102dc6bd`.

| Quantity | Approximate directed result |
|---|---:|
| Reused mixed bound's lambda-l1 upper |18.61447623568004556|
| Signed kernel cross nu0 |+0.181002566603288306|
| Sharp forward error upper |5.304632297062084e−69|
| Physical Gram cross error upper |1.669239650964427e−66|
| Complete new transfer allowance |4.112478850097329e−11|
| Actual theta−mu interval |[0.18100256656216351,0.18100256664441310]|
| Actual theta interval |[−4.331313289332710,−4.325974508801690]|
| Frozen Clin+2theta interval |[16.865709521960626,16.876387083066322]|

Displayed intervals round outwards; the saved integer endpoints and
exact rational allowances control all tests. Reopening a saved interval
as an Arb ball conservatively inflates its width. This includes the
reused proxy interval and does not change its source or error contract.
The actual correction is positive under the inherited certificates,
while the unchanged sufficient test needs a negative upper below about
−8.25719. Even its zero-transfer-error expression remains above16.86.
Further reduction of this transfer allowance cannot make the fixed
sufficient expression negative. This determines a limitation of the
frozen envelope, not the sign of the actual center below that envelope.

The useful next obligation is to derive a tighter linear upper or
restore discarded signed contributions, preserving their common
realization and the old certificate; replaying this mixed correction
is unnecessary.
The good-block metric, actual full-border sign, retained/cofinal
positivity, Robin and RH remain unresolved. All inherited actual-form,
domain, source/lift and component premises remain conditions. This is
a numerical application of existing bounds, with no new generic
algorithm, literature novelty or Lean theorem claim.

For an explicitly requested reproduction only, choose a fresh output:

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_signed_correction.py --out /tmp/c27-signed-correction-reproduction.json.zlib
```

## All-real steering capacity from saved scalar certificates

[full_gamma_steering_capacity.py](full_gamma_steering_capacity.py) consumes
only the saved negative-channel, coherent-mixed and signed-correction
scalar packets. The source, core, lift, far row, legal t0 and frozen
fee D* remain fixed. It bounds the comparison expression
E*(alpha)=alpha²/t0+D*+(alpha*rho+betaJ+f)Qtrue(alpha*rho+betaJ+f)*
for every real alpha. This algebraic domain does not assert that every
alpha is a legal physical test; any legal subset with the same fixed
data is included in the exclusion.

The existing identity rho0²/delta=1/t0+rr−H11 and Cauchy give a uniform
absolute derivative cap, using the saved far positive energy A and
negative-coordinate radius. The existing lower inverse comparison
pays eta/(1+eta) times the positive reference norm. A two-step triangle
bound encloses that norm by4*r0+4*s²*rr+2*A, with s=alpha−alpha0.
Scalar completion then gives E*(alpha)>=K−B²/H for all real alpha.
The finite lower is read directly from its saved scalar quadratic;
subtracting upper allowances from Clin would not justify that lower.
The complete proof and common-object conditions are in the owning
[application exposition](https://github.com/the-omega-institute/trureturing/blob/dev/Library/Weil/connesconsanimoscovici2026spectral.md#test-all-steering-coefficients-before-a-new-signed-supplier).

The sole fixed1024/640 application takes0.61818seconds including
serialization and uses60 unchanged source/input bindings. No old
producer, new derivative, source, norm, feature, matrix, inverse or
tail computation is run.

| Quantity | Approximate directed value |
|---|---:|
| Negative-square credit upper |2.327429771347697837|
| Absolute signed-derivative capacity upper |515.620424366832531834|
| Comparison curvature lower |25136.371399957733363764|
| Comparison base lower |10.639254031103735470|
| Maximal steering gain upper |10.576881515391657423|
| All-alpha comparison floor lower |0.062372515712078048|

The saved integer lower endpoint exceeds0.062372515712078048.
Consequently this fixed-fee/t0 comparison recipe cannot become
negative for any alpha; an actual derivative supplier for this recipe
is unnecessary. D* is a frozen charge upper used to define this
comparison family. Its positive floor is not a lower bound on an
actual center with a smaller fee. At alpha0 the separate self-term-only
refinement class retains a floor above10.639254031103735470; it bounds
that frozen upper-expression class, not the actual center.

The [directed packet](full-gamma-steering-capacity.json.zlib) has SHA256
`85a52adf97e96def95aaf5abaa307338522a2bae5d7b40c96692ad0366a444fe`;
the source hash is
`dcdb47a55003e4ca745bf33ce81133fd65577facdf6004d38308fb9f1dddc96e`.
This is a source-specific application of existing Schur, Cauchy and
quadratic estimates, with no new generic theorem, algorithm, literature
novelty or Lean claim. All inherited actual-form/domain, common-source,
short/lift, weight and component premises remain. Actual-center and
full-border signs, good-block metric, retained/cofinal positivity,
Robin tail and RH remain unresolved. A next route must change a
justified fee, coefficient, source or comparison rather than supply
the derivative of the excluded recipe.

For an explicitly requested reproduction only, choose a fresh output:

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_steering_capacity.py --out /tmp/c27-steering-capacity-reproduction.json.zlib
```

## Mixed-H norm and linked fee/eligibility repricing

[full_gamma_mixed_norm_reprice.py](full_gamma_mixed_norm_reprice.py)
consumes the signed-correction packet's saved112 left features a=P0*ell
and the original physical112-by112 kernel K0. One new quadratic bounds
the weighted mixed energy by a^T*K0*a+epsfar*Zlambda². The original
physical epsfar retains its512-bit metadata units and K0 already
contains the far weight; no extra kappa or weight enters this norm.
The kernel is not assumed PSD. No source, feature, moment, Gram,
inverse or signed-response supplier is repeated.

The same total energy upper min(Eold,(P+epssm+mixedUpper)²) enters both
the independent source-valid fee Kconst_up+kappa*E and the complete
reflected Bessel budget. The resulting legal coefficient is
tnew=[hlo/(hlo+p_up)]/hhi. The shifted reference curvature is
H11new=H11old+1/tnew−1/told; H12, rho, betaJ, f, r0, rr and theta0
remain fixed. Positive/negative Cauchy caps scale with E and sqrt(E),
respectively, using their original homogeneous formulas. The existing
steering capacity function is reused at these new common parameters.
The complete derivation and conditions are in the owning
[application exposition](https://github.com/the-omega-institute/trureturing/blob/dev/Library/Weil/connesconsanimoscovici2026spectral.md#reprice-the-same-source-from-its-saved-mixed-h-features).

The sole fixed1024/640 application takes0.55967seconds including
serialization, uses62 unchanged runtime bindings and makes one new
mixed quadratic. It reuses all old science and accepted tests.

| Quantity | Previous value | Repriced value |
|---|---:|---:|
| Weighted mixed-H norm upper |0.478756042961117000|0.052637674256435558|
| Weighted far-energy upper |9.164839739390695|6.766398882441490|
| Positive source-valid fee |178.644404880862231|176.171012747133362|
| Sufficient eligibility cap |0.035871186206647585|0.042365798558350576|
| Legal coefficient |0.000033349857564549898|0.000039387973940699242|
| Fixed-point sufficient expression |[16.865709521960626,16.876387083066322]|[6.738778327421815964,6.749455888527510587]|

Norm/energy/fee displays round upward and legal cap/coefficient displays
downward; intervals round outward. Exact saved grid endpoints and
rationals are authoritative. The new finite derivative is approximately
−174.6152875359602151, and the all-alpha capacity lower floor is near
−16.43212445437972714. This new recipe is not excluded by that capacity
test. No favorable actual derivative or negative center follows: the
fixed-point sufficient upper remains positive. Its separate self-term
refinement class retains a positive floor above2.141787799175894749;
that class diagnostic is not an actual-center lower bound.

The old fixed-fee/t0 family exclusion remains true. This package
supplies one improved source-energy enclosure and linked comparison,
not an optimized all-energy family or a new generic theorem/algorithm,
literature-novelty result or Lean proof. All actual-form/domain,
common-source/short/lift, original-weight, homogeneous-cap and component
premises remain. Actual-center/full-border signs, good-block payment,
retained/cofinal positivity, Robin tail and RH remain unresolved.
The next steering datum is f*Qref*rho*, or a sufficient conditional
capacity enclosure using the known theta0 and saved full physical Gram;
neither has been supplied by this packet.

The [directed packet](full-gamma-mixed-norm-reprice.json.zlib) has SHA256
`66360516946c571c187b5a4eb35aeed1487eef9eff28bee26d5f8fb3211f1a59`.
Source SHA256 is
`92d38517ef0e724cc819ae26e934f492273d93a2a8634c89b573f98c3154b518`.
For an explicitly requested reproduction only, choose a fresh output:

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_mixed_norm_reprice.py --out /tmp/c27-mixed-norm-reprice-reproduction.json.zlib
```

## Saved signed-response joint bound and energy-only family

[full_gamma_negative_index_joint.py](full_gamma_negative_index_joint.py)
consumes saved scalar certificates only. It reuses the classical
negative-index-one Schur inequality and scalar square completion; no
source, feature, norm, signed pairing, derivative, Gram, inverse, trace
or tail supplier is evaluated.

Use the SAME $P=TRT^*\succeq0$, $Q_R=P-e_0e_0^*/\delta$,
$y_0=\alpha_0\rho+\beta_J$, far row $f$ and signed
$\vartheta_0=fQ_Ry_0^*$. For saved $\eta_{\rm up}<1$, set
$\varepsilon=\eta/(1+\eta)$ and
$Q_\varepsilon=(1-\varepsilon)P-e_0e_0^*/\delta\preceq Q_{\rm true}$.
The saved whole interval $q_0$ and old finite derivative $d_0$ give

$$
\begin{aligned}
b&=D_{\rm old}+\alpha_0^2/t_{\rm old}-q_0+\varepsilon r_0,\\
g&=d_0-\alpha_0/t_{\rm old}-\varepsilon\rho Py_0^*,\\
c&=q_0-D_{\rm old}+\vartheta_0-\alpha_0d_0
-\varepsilon(\beta_J+f)Py_0^*.
\end{aligned}
$$

Let $0<B\le b$, $G\ge|g|$, $M\ge|c|$ use the Cauchy bounds in the
source. The positive far-row energy UPPER in these constants is frozen
at the accepted new cap, including at the artificial energy LOWER
endpoint. For every real $\alpha$, the paid comparison is bounded below
by $D-kM^2/(Bk-G^2)$ when $Bk-G^2>0$, where $k=1/t$.

Two specified interfaces are evaluated once: the accepted new fixed
fee/coefficient, and an endpoint for the entire source-energy-only
pricing rule. The latter uses the saved signed Hankel cross to obtain

$$
E_H=\left[\sqrt{\max(w_F\operatorname{cross}_{\rm lo}/4,0)}
-\epsilon_{\rm sm,up}-H_{\rm new,up}\right]_+^2\le E_{\rm true}.
$$

The omitted exact principal diagonal is nonnegative. All three terms
must use the same weighted norm. With the charge constant, original
Bessel coefficients and $h$ endpoints fixed, both $D(E)$ and $k(E)$
are nondecreasing in the source energy upper $E\ge E_{\rm true}$.
An unrounded-budget reciprocal coefficient at $E_H$, rounded downward,
therefore provides a lower endpoint. A positive endpoint floor excludes
all real-$\alpha$ comparisons obtained solely by reducing that energy
upper in this fixed pricing rule. It does not cover better independent
charge terms, different budget formulas, sources or lifts. The formal
endpoint need not itself be a legal source-energy cap.

The fee is an upper defining the paid comparison, never a lower bound
for an actual smaller charge. All inherited actual-form/domain,
source/short/lift, reference inverse, component and physical-weight
premises remain. Actual-center/full-border signs, good-block payment,
retained/cofinal positivity, Robin and RH remain unresolved. There is no
new generic theorem, algorithm, novelty or Lean-proof claim.

The single fixed1024-working/640-export application uses no new matrix,
norm, source, feature, derivative, signed-response, inverse, trace or
tail supplier. Its64 source/input bindings agree at entry and prewrite;
18 scalar intervals and14 exact rationals are retained. Stored assembly
time is0.15059seconds; total time including serialization is0.18470seconds.

| Bound or comparison | Directed value |
|---|---:|
| Timelike lower $B$ |>200.284690047319663558|
| Pairing upper $G$ |<1139.054295521911259039|
| Pairing upper $M$ |<161.346407789081724764|
| Source-energy lower $E_H$ |>6.227486948459329823|
| Fixed repriced rank denominator |>3787475.149285169293|
| Fixed repriced all-real floor |>1.667100663167811801|
| Energy-endpoint rank denominator |>3579664.592067604194|
| Energy-endpoint all-real floor |>-1.473501281813134057|

The current repriced fee/coefficient comparison is therefore strictly
positive for EVERY real steering coefficient under its inherited
interfaces. Computing a favorable derivative cannot make that specified
paid comparison negative. This strengthens the preceding independent-cap
test, whose negative lower floor left the same family unexcluded.
It still supplies no lower bound for an actual smaller-fee center.

The energy endpoint has a positive rank denominator, but its saved floor
is negative. The sufficient endpoint test therefore does NOT exclude the
whole fixed-budget energy-cap-only family. That negative lower bound is
neither a realizable negative value nor a negative upper. The accepted
actual energy now lies between the certified lower6.227486948459329823
and upper6.766398882441490; a sharper SAME-source energy enclosure or a
change to another charge/budget component can escape the scope of the
fixed-recipe exclusion. No such escape or successful negative recipe is
certified here. In particular independent triangle bounds do not supply
the signed principal/mixed correlation needed to sharpen that enclosure.

The [source](full_gamma_negative_index_joint.py) and [exact packet](full-gamma-negative-index-joint.json.zlib)
retain packet SHA256
`a7df6338b629e9dd791ac2914c8efc7f3ef3c69c20e562da0aa18067aa40728d`
and source SHA256
`9c84ee8891ccc7d9af9aeaa79b162e3c0854bf714badd9366607890aa8e345c3`.
The exact packet endpoints and rational comparisons are authoritative.
The two scalar applications use one common realization and frozen far-row
upper; they are not a parameter scan. The actual good-block metric,
full-border/cofinal positivity, Robin and RH remain unresolved. No Lean
verification or new generic theorem is claimed.

Reproduction in the pinned environment, with a fresh output:

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_negative_index_joint.py --out /tmp/c27-negative-index-joint-reproduction.json.zlib
```

## Shared principal/mixed projections in the source energy

[full_gamma_shared_projection_energy.py](full_gamma_shared_projection_energy.py)
reuses the saved signed principal and mixed projections on the SAME
physical representer b=Cphys t. One112-component quadratic h^T K0h,
h=P0t, and the original physical operator error give

$$
\|b\|_{W_F}^2\le h^TK_0h+\epsilon_{\rm far}K_{\rm phys}L_V^2.
$$

Take the smaller outward upper with the previous direction-norm cap
squared. The kernel is already weighted; no extra far weight or kappa
is inserted, and K0 need not be PSD. Kernel arrays use640-grid units;
original physical error and Kphys metadata use512-grid units.

In uF=v+r+c, the saved mu/kappa plus its separate sharp forward
allowance encloses s=<v,b>; the saved mixed signed cross/kappa plus its
physical allowance encloses z=<c,b>. The smooth remainder r is not in
z. The classical three-vector Gram/Schur inequality bounds <v,c> by
s*z/N2 plus or minus the product of the two residual-norm caps for any
N2>=||b||². Only a norm UPPER is needed. Interval evaluation and ordinary
Cauchy supply joint bounds Jlo<=||v+c||²<=Jup, uniformly for the same
realization, without asserting simultaneous attainability of endpoints.

The accepted pricing maps then consume the valid regrouping (v+c)+r+0.
The final zero is a component of this new decomposition, not a claim
that the original c vanishes. The energy lower is the reverse-triangle
bound [sqrt(Jlo)-smooth]+², strengthened by the old lower if necessary.
The smaller energy upper feeds BOTH the fee and the full original
Bessel budget. The saved B/G/M and their fixed far-row UPPER remain
unchanged at the fixed and energy-lower comparisons.

The single fixed1024-working/640-export application evaluates one
physical representer quadratic. Its66 runtime bindings agree at entry
and prewrite. Stored assembly time is0.61203seconds; time including
serialization is0.72872seconds. No old source, feature, mixed norm,
Gram, inverse, signed projection, trace or tail producer is replayed;
there is no source/parameter/precision/cutoff scan.

| Quantity | Directed value |
|---|---:|
| Representer norm upper |<3.918738184778600392|
| Principal/mixed cross interval |[−0.113379023101042651,0.013416593112150290]|
| Source-energy lower |>6.268984281541914933|
| Source-energy upper |<6.524928095642312689|
| Source-valid positive fee |<175.921995998246711018|
| Legal coefficient |>0.000040122278284528281535|
| Fixed repriced all-real floor |>0.296481006612169882|
| Energy-endpoint all-real floor |>−1.220854294197841680|

The previous energy bracket was [6.227486948459329823,
6.766398882441490]. Both ends improve. The new fixed fee/coefficient
comparison nevertheless stays positive for every real steering
coefficient. Further derivative supply cannot make that specified paid
recipe negative. This is not an actual smaller-fee center lower bound.
The negative energy-endpoint lower remains inconclusive for the whole
fixed-budget energy-cap-only family; it is neither a negative upper nor
a realized value. The cross interval still contains positive numbers,
so this packet does not certify a negative principal/mixed cross.

The [exact packet](full-gamma-shared-projection-energy.json.zlib) has
SHA256 `59374362635382b3b180aff21584a5a1e32d603b19ccaa5db3ea693086b5f89c`.
Source SHA256 is
`ba5ea56e4a30c3afcaf52d428bbd060f19ca18a29464f6dc15e8055837329735`.
These are applications of existing Gram/Schur and physical-error
certificates, not a new generic theorem, algorithm, novelty or Lean
proof. Actual-form/domain, source/short/lift, original-weight, component
and inverse-comparison premises remain. Actual-center/full-border signs,
good-block payment, retained/cofinal positivity, Robin and RH remain open.

Reproduction only when explicitly needed, using a fresh output:

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_shared_projection_energy.py --out /tmp/c27-shared-projection-energy-reproduction.json.zlib
```

## Directional full-jet norm from three retained moments

[full_gamma_directional_jet_norm.py](full_gamma_directional_jet_norm.py)
pays omitted moments on the SAME saved real feature vector h=P0t.
It reuses the saved Q=h^T K0h and the sharp weighted physical-to-jet
forward allowance; no quadratic or physical supplier is repeated.

For each omitted n=k+l>=56, it groups the complete ordered signed
feature sums a_n=sum A_k A_l, g_n=sum A_k B_l and b_n=sum B_k B_l
before bounding a_n U_(2n+2)-2g_n T_(2n+3)+b_n V_(2n+4). The first
upper uses the positive parts of a_n,b_n and absolute g_n with the
original moment integral caps. Ordered AA/BB terms and the one AB/BA
factor2 are retained.

The second upper reads only the already retained U112,T113,V114
intervals. Their TRUE two-by-two PSD Gram Sstar dominates every omitted
moment block Sn by Sn<=q^(2(n-55))*Sstar, q=32768/32769. The exact
positive-part trace support formula for an indefinite coefficient matrix
bounds its contraction. Its interval implementation uses tauUP and
deltaLOW, not a PSD assertion about rounded centers. The per-n minimum
of both valid uppers is uniform on the SAME realization.

The accepted outer certificate separately pays the retained moment
midpoint error112*2^-640 and alpha midpoint error using the original
explicit TRACE upper Hstar. Entry masking decreases Frobenius norm;
PSD of the full true Gram then bounds the masked norm by Hstar.
The masked midpoint or omitted block need not be PSD. The original
512-grid complete physical error is not relabeled as a640-grid moment
error. Adding the old sharp forward norm once gives the full physical
representer cap without an extra far weight or kappa.

A strictly improved cap alone enters the accepted projection/pricing
maps. The valid regrouping is still (v+c)+r+0. Stronger old energy
endpoints and unchanged comparison outputs are reused. The SAME saved
B/G/M and fixed far-row energy UPPER are used throughout.

The single fixed1024-working/640-export application retains68 runtime
bindings, unchanged at entry and prewrite. The saved quadratic is reused.
The omitted-block payment is below$1.43015\cdot10^{-130}$ and the
separate rounding payment below$1.25068\cdot10^{-185}$. Adding the
original sharp forward allowance gives the following outward bounds.

| Quantity | Directed value |
|---|---:|
| Representer squared-norm upper |<11.172943110390517643|
| Representer norm upper |<3.342595265716523695|
| Principal/mixed cross interval |[−0.076771255960990555,−0.060620829249211873]|
| Source-energy lower |>6.342199815821979165|
| Source-energy upper |<6.376853250919508462|
| Source-valid positive fee |<175.769293814626319160|
| Legal coefficient |>0.000040586603433261736126|
| Fixed repriced all-real floor |>−0.572919445838071768|
| Energy-endpoint all-real floor |>−0.779703303246544574|

The norm cap improves from3.918738184778600392, and both ends of the
preceding energy bracket [6.268984281541914933,6.524928095642312689]
tighten. The principal/mixed cross is now strictly negative under the
retained common-source premises. This follows from the sharper common
direction norm, with no new signed pairing supplier.

The preceding paid-recipe exclusion applies to its original fee and
coefficient. The newly repriced pair has a negative all-real LOWER,
so that lower-bound test leaves the new recipe unexcluded. Its negative
lower neither realizes a negative comparison nor supplies a negative
upper. The energy-endpoint test also remains inconclusive for the whole
fixed-budget energy-cap-only family. Actual smaller-fee center sign is
not determined. Displayed endpoints round outwards; exact rationals
and640-bit packet endpoints are authoritative.

The [exact packet](full-gamma-directional-jet-norm.json.zlib) has SHA256
`9ccd6494c0520f9f53efe0b44bcae879566157c7e7ffb698165f9d02c258fd41`.
Source SHA256 is
`d9cbc7f94640ec60a885dfce018fc626726fe77caaaa87aef63c912c0f70abcd`.
Stored assembly time is0.14732seconds; time including
serialization is0.18160seconds. There is one new grouped omitted-block
contraction, three signed feature-antidiagonal families, three anchor
scalar reads and one feature absolute norm. One projection, one pricing
and two rank-floor applications consume the new cap; all old producers,
quadratics, signed pairings, moments, sources, lifts, inverses, traces
and long-tail suppliers are reused. There is no source, parameter,
precision or cutoff scan.

This applies existing jet, moment, semidefinite support-function and
midpoint estimates. No new generic theorem, algorithm, novelty or Lean
proof is claimed. The actual-form/domain, common-source/short/lift,
original-weight, component and inverse-comparison premises remain.
Actual-center/full-border signs, good-block payment, retained/cofinal
positivity, Robin and RH remain unresolved.

Reproduction only when explicitly needed, using a fresh output:

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_directional_jet_norm.py --out /tmp/c27-directional-jet-norm-reproduction.json.zlib
```

## Equality row inside the signed-response and energy relaxation

A negative lower bound does not establish that its relaxation permits a
negative value. Test one explicitly constructed row before replacing the
accepted timelike estimate by a larger constrained optimization. This is
an application of classical timelike Schur equality on the SAME finite
reference form, with no new generic inequality or optimality claim.

Retain $P=TRT^*\succeq0$, $Q_R=P-e_0e_0^*/\delta$, and
$y_0=\alpha_0\rho+\beta_J$. The existing scalar packet already encloses

$$
r_0=y_0Py_0^*,\quad r_\rho=\rho P\rho^*,\quad
r_{\rho\beta}=\rho P\beta_J^*,\quad
\lambda=y_0[0],\quad g=\rho Q_Ry_0^*.
$$

Use the saved finite cross, rather than supplying another Gram:

$$
p=\alpha_0r_\rho+r_{\rho\beta}=\rho Py_0^*,\qquad
b=\lambda^2/\delta-r_0=-y_0Q_Ry_0^*>0.
$$

Both the coefficient $\lambda$ and the finite parameters retain their
interval uncertainty. No rounded center is declared the true form.
Choose the one rational steering value $\alpha_*=101/2000$, and an
exact interior value $\theta_*$ of the saved signed-response interval.
The chosen value is not asserted to equal the actual signed response.
Define, for the true common finite parameters,

$$
s=\alpha_*-\alpha_0,\qquad
\tau=-\frac{\theta_*+sg}{b},\qquad
f_{\rm rel}=\tau y_0-s\rho.
$$

Then, exactly,

$$
f_{\rm rel}Q_Ry_0^*=\theta_*,\qquad
\alpha_*\rho+\beta_J+f_{\rm rel}=(1+\tau)y_0.
$$

This is a row in a selected scalar relaxation. Its generation by the
original physical source/coupling, and membership of $\alpha_*$ in a
physical legal subset, are not supplied. The construction defines a row
for each actual finite reference form covered by the common intervals.

Its two norm/cap checks use

$$
f_{\rm rel}Pf_{\rm rel}^*=
\tau^2r_0-2\tau sp+s^2r_\rho,\qquad
f_{\rm rel}[0]=\tau\lambda-s\delta\frac{p-g}{\lambda}.
$$

The second identity uses $\lambda\ne0$ and
$\rho[0]=\delta(p-g)/\lambda$. The original far-row positive-energy
and negative-coordinate bounds are homogeneous in the source-energy
UPPER. From the saved mixed-norm recipe's $E_0,A_0,C_0$ and the latest
source-energy upper $E$, the accepted interface gives outward caps

$$
A=A_0E/E_0,\qquad C=C_0\sqrt{E/E_0}.
$$

These are bounds in different units: $E$ is source energy, $A$ bounds
$fPf^*$, and $C$ bounds $|f[0]|$. No source-energy lower is substituted
for a far-row upper. Check the constructed row against BOTH caps.

The already certified upper inverse comparison gives
$Q_{\rm true}\preceq Q_R+\eta(1-\eta)^{-1}P$. Thus this constructed
row's paid comparison has the sufficient UPPER

$$
\mathcal E_{D,k}^{\rm rel}(\alpha_*)
\le k\alpha_*^2+D
+(1+\tau)^2\left(-b+\frac{\eta}{1-\eta}r_0\right).
$$

It is an upper on this particular constructed row, rather than the
previous all-real lower bound. All scalar intervals propagate, and no
matrix, inverse, physical far row or actual derivative is regenerated.

The single application satisfies both caps and has a strictly negative
UPPER for this constructed row:

| Quantity | Directed value |
|---|---:|
| Constructed positive-energy upper |<0.583108075508182745|
| Allowed positive-energy cap |>2.712924125459587053|
| Constructed absolute first-coordinate upper |<2.228217728848581443e-16|
| Allowed absolute first-coordinate cap |>4.367163575891468984e-16|
| Constructed paid-comparison upper |<−0.567072019035610609|

All three checks, and their conjunction, pass. Consequently these
selected scalar constraints do permit a negative paid comparison under
the retained finite-form and inverse-comparison premises. No uniform
nonnegative bound follows from this relaxation alone. This is a
feasible row in the relaxation, not just a negative lower bound and not
an attained value of the original physical source. The chosen response
lies inside the measured interval but is not identified with its true
value. Displayed values round outwards; exact rational caps and640-bit
interval endpoints are authoritative.

The [source](full_gamma_relaxation_row_witness.py) and
[exact packet](full-gamma-relaxation-row-witness.json.zlib)
retain packet SHA256 `1879f4b0653e762ca0d7bb1e1933a59afc79fb10851d155f1c16c87055872992` and source SHA256 `218946086714e3408a4837a98e211622a10e7361eb0789427b1d3c9f24245340`.
The single fixed1024-working/640-export application uses one literal
$\alpha_*$ and one interior target, with no parameter/source/precision/
cutoff scan or retry. Earlier producers and passed fixtures are reused.

The constructed row lies in the span of $\rho,y_0$, so a new payment
confined to the perpendicular residual would vanish on this row.

The selected relaxation keeps the finite form, inverse comparison,
signed-response interval and two homogeneous far-row caps. It omits the
full original source-to-row image constraint and other physical
observations. A row admitted by this relaxation is not automatically
admitted by the original source model. Determining the actual-source
sign requires its response or an additional relation that distinguishes
it from the constructed row; this test does not establish which single
additional scalar would suffice.

This is a numerical application of existing equality and inverse
comparison identities, with no new generic theorem, novelty or Lean
proof. Actual-center/full-border signs, good-block payment,
retained/cofinal positivity, Robin tail and RH remain unresolved.


The saved assembly time is<0.32502seconds; time including serialization
is0.42767seconds. Runtime guards bind11 invoked sources and59 inherited
sources/data, unchanged at entry and prewrite. One constructed-row
application and one homogeneous cap rescaling are new. No actual-source,
derivative, Gram, inverse or tail supplier is evaluated, and no earlier
producer or accepted fixture is replayed.

Reproduction only when explicitly needed, using a fresh output:

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_relaxation_row_witness.py --out /tmp/c27-relaxation-row-reproduction.json.zlib
```

### Retain the original source on a zero-first-coordinate direction

The scalar relaxation omits the full source-to-row relation. Reuse
$P=TRT^*\succeq0$, $Q_R=P-e_0e_0^*/\delta$,
$y_0=\alpha_0\rho+\beta_J$ and $\lambda=y_0[0]\ne0$, and set

$$
c=\rho[0]/\lambda,\quad w=\rho-cy_0,\quad w[0]=0.
$$

Then $Q_Rw^*=Pw^*$ and
$fQ_R\rho^*=\zeta+c\theta$, with
$\zeta=fPw^*$ and the SAME actual signed response
$\theta=fQ_Ry_0^*$. The original full row is
$f=\kappa u_F^*W_F\mathsf B$, $\mathsf B=R_FV$,
$\kappa=33/32$. Thus for $b=\mathsf BQ_Ry_0^*$ and
$d=\mathsf BPw^*$,

$$
\theta=\kappa\langle u_F,b\rangle_{W_F},\qquad
\zeta=\kappa\langle u_F,d\rangle_{W_F}.
$$

Keep the whole actual-$\theta$ interval, the source-energy UPPER
$E\ge\|u_F\|_{W_F}^2$, the saved sharper $J_b$, and the common
physical Gram relation. Source energy $E$ differs from the far-row
positive-energy cap $A$. The existing `projected_energy` bound is
applied to the SAME $(u_F,d,b)$ with caps $(E,J_d,J_b)$; its output
is intersected with the old scalar-only bound
$|\zeta|\le\sqrt{AwPw^*}$. Empty intersections stop.

The [coarse consumer](full_gamma_zero_channel_pairing.py) performs
one cached reference-vector action, one112-by18 feature action, one
new feature quadratic and one old/new feature cross. Its uniform
physical error pays over2984.69 against a quadratic of13.75086838,
so it supplies neither strict narrowing nor constructed-row exclusion.
The old direction quadratic and all original source/sine/moment/tail/
inverse/solve/trace producers are reused without replay.

The [directional consumer](full_gamma_zero_channel_jet_refine.py)
reuses those saved new features, quadratic and cross with the existing
`directional_norm` method. In the original sharp forward allowance
$\epsilon_b^\sharp\ge C_{\rm forw}L_b$, all coefficient factors are
direction-independent for the SAME physical operator, support2048,
scale32768, order56, symbol, moment, weight and rounding certificates.
Since $L_b>0$, outward rounding gives

$$
\epsilon_d^\sharp=\operatorname{up}
 (\epsilon_b^\sharp L_d/L_b)\ge C_{\rm forw}L_d.
$$

This transfers a uniform error coefficient, not an arbitrary actual
norm cap. The original three retained anchors and existing omitted-jet
and midpoint payments bound $J_{h_d}$; select the smaller of the coarse
cap and $\operatorname{up}((\sqrt{J_{h_d}}+\epsilon_d^\sharp)^2)$.
Apply the SAME common-projection bound with the saved old/new cross,
then intersect with the coarse interval. No direction, feature,
quadratic, cross or moment supplier is regenerated.

| Quantity | Outward displayed enclosure or cap |
|---|---:|
| Coarse squared-norm upper |<2998.447008|
| Directional squared-norm upper |<13.750868382|
| Transferred forward-error upper |<1.416877e-67|
| Retained scalar-only $\zeta$ interval |contained in[−66.425463,66.425463]|
| Refined actual $\zeta$ interval |contained in[−9.656778,9.656778]|
| Constructed row's $f_{\rm rel}Pw^*$ |contained in[−20.900041,−20.900040]|
| Refined actual $fQ_R\rho^*$ interval |contained in[−33.775533,−14.432247]|

The norm strictly improves and the actual zero-pairing interval
strictly narrows; it excludes the particular saved constructed row
under the retained original-source premises. It does not exclude all
negative rows or supply the sign of $\zeta$. The reference-$\rho$
pairing is strictly negative in this conditional model, but a far-row
pairing is neither an exact derivative nor the paid expression.
Actual-center/full-border signs, good-block payment, retained/cofinal
positivity, Robin tail and RH remain unresolved.

The [coarse packet](full-gamma-zero-channel-pairing.json.zlib)
has SHA256 `c4cbbead777f652812c0aa0bd238bf08053daa52be7afbedfa1321da19847a9b`;
its source SHA256 is `612ba841171f1a040c90f376c85cbc326f3e46c2a9ef94e0257094720c0460bb`.
The [refinement packet](full-gamma-zero-channel-jet-refine.json.zlib)
has SHA256 `ee401d13ef503a44089c33e958b8f90dab09aefa843050f09713ee6a33729861`;
its source SHA256 is `4c07a0a9a55dd9d87c1f1b741960f29c495246a9a868ecd3b2341bf0b11c07d0`.
Exact rational caps and640-bit endpoints are authoritative; the
displayed actual-pairing bounds also enclose the exported Arb balls.

Each consumer has one fixed1024-working/640-export application, with
no parameter/source/precision/cutoff scan or retry. The refinement
adds one homogeneous-forward, one directional-norm and one projection
application on saved output, with zero repeated finite actions.
Coarse entry/prewrite guards bind19 invoked and53 inherited paths;
refinement guards bind19 invoked and55 inherited paths. Saved assembly
times are<0.472150seconds and<0.703493seconds respectively; these
exclude serialization. No generic theorem, novelty or Lean claim is
made. Original actual-form/domain, source/lift and physical-error
certificates remain premises.

Reproduction only when explicitly needed, using fresh outputs:

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_zero_channel_pairing.py --out /tmp/c27-zero-channel-reproduction.json.zlib
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_zero_channel_jet_refine.py --out /tmp/c27-zero-channel-refinement-reproduction.json.zlib
```

### Signed old/new cross from the saved directional features

The [bilinear consumer](full_gamma_zero_channel_cross_refine.py) uses
the SAME actual-symbol full jets, saved features and masked cross from
the preceding two consumers. It supplies a two-sided physical-cross
payment without repeating any finite direction, feature, quadratic,
cross or original producer. Approximate-symbol finite proxy samples
are not substituted for these full jets.

For each omitted antidiagonal $56\le n=k+l\le110$, sum ordered AA,
AB+BA and BB feature products into $a_n,g_n,b_n$. The two-sided
omitted payment is bounded upward by

$$
w_F\alpha_{\rm up}\sum_{n=56}^{110}
\left(\frac{S^2|a_n|}{2n+1}
+\frac{S|g_n|}{2n+2}+\frac{|b_n|}{2n+3}\right).
$$

Absolute AA and BB coefficients pay both signs; AB+BA already includes
both mixed orientations. Feature uncertainty is included before
taking absolute uppers. The original entrywise moment intervals and
full-moment trace cap pay the masked midpoint/prefactor error by

$$
w_F\sqrt{H_bH_d}
\left[\alpha_{\rm up}e_H+r_\alpha(H_*+e_H)\right].
$$

The original physical-vector bounds pay a further

$$
\epsilon_b^\sharp\sqrt{J_{j_d}}
+\epsilon_d^\sharp\sqrt{J_{j_b}}
+\epsilon_b^\sharp\epsilon_d^\sharp.
$$

Add these payments to the saved masked-cross interval and intersect
with its accepted physical interval. No masked or rounded matrix is
assumed PSD. The existing common projection uses the new signed
cross, the entire actual $\theta/\kappa$ interval and the SAME
source/physical norm bounds. Retain intersections with both preceding
pairing intervals; empty intersections stop the application.

| Quantity | Outward enclosure or directed bound |
|---|---:|
| Bilinear masking payment |<1.055086e-128|
| Bilinear moment/prefactor rounding payment |<2.575728e-185|
| Total additional bilinear cross-error payment |<4.932753e-67|
| Actual physical old/new cross |contained in[10.6404166919614235901773440,10.6404166919614235901773441]|
| Actual $\zeta$ |contained in[−8.423007,0.178344]|
| Actual $fQ_R\rho^*$ |contained in[−32.541761,−23.910682]|
| New fixed-recipe all-real lower, rounded down |−2.756762331|
| New energy-endpoint all-real lower, rounded down |−2.959238524|
| Retained stronger fixed-recipe lower, rounded down |−0.572919446|
| Retained stronger energy-endpoint lower, rounded down |−0.779703304|

The displayed pairing bounds also enclose the exported Arb balls.
The final physical-cross interval retains the saved masked-cross width
of about $5.56\times10^{-27}$; the additional payment is not its final
width. Physical-cross and source-pairing strict-narrowing flags are
true. The $\zeta$ interval contains zero. The negative reference-$\rho$
pairing does not supply a paid-expression sign.

The existing all-real square completion is applied with explicit
translations $k-k_0,D-D_0$ in the old reference base, curvature and
slope. It retains the SAME-row old negative-square credit and
inverse-loss uppers. The fixed current recipe and saved energy-lower-
endpoint recipe are tested separately. Both new exclusion flags are
false; their lowers are weaker than the preceding timelike joint
bounds, which remain the stronger certificates. A nonpositive lower
constructs no negative row or actual-source upper. A positive endpoint
lower would exclude only the fixed monotone energy-cap-only pricing
family, not an actual center with a smaller justified charge.

The [exact packet](full-gamma-zero-channel-cross-refine.json.zlib)
has SHA256 `5cfba16d5428c1c61618a861a8167ca6828f46142d7bd7993d4d84b20de6e980`;
source SHA256 is
`a35c3222a788508622b92adc6525bc73a26cdc2070c26eb5a701a62015e991a5`.
One fixed1024-working/640-export application performs one bilinear
payment, one existing common projection and two scalar floor
applications. There are18 invoked and58 inherited bindings, with zero
repeated finite actions or original producers and no parameter/source/
precision/cutoff scan or retry. Saved assembly time is<0.074813seconds,
excluding serialization. Original actual-form/domain, common-source/
lift and physical error premises remain. Actual-center/full-border/
good-block/cofinal positivity, Robin and RH remain unresolved; no
generic theorem, novelty or Lean proof is claimed.

Reproduction only when explicitly needed, using a fresh output:

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_zero_channel_cross_refine.py --out /tmp/c27-zero-channel-cross-reproduction.json.zlib
```

### Both source responses in the paid indefinite plane

The [joint-plane consumer](full_gamma_joint_plane_paid_floor.py) reuses
the saved finite scalars and both signed responses. The classical
two-dimensional Schur projection retains a positive square that the
single-timelike lower discards. No finite or physical supplier is
repeated, and the original actual-symbol certificates remain premises.

On the SAME real form $Q_\varepsilon=\gamma P-e_0e_0^*/\delta
\preceq Q_{\rm true}$, use

$$
b=-y_0Q_\varepsilon y_0^*>0,\quad
g=\rho Q_\varepsilon y_0^*,\quad h=\rho Q_\varepsilon\rho^*,
\quad t=wQ_\varepsilon w^*\ge0,
$$

with $w=\rho-cy_0$, $w[0]=0$, $\omega=c+g/b$ and
$a=t+b\omega^2=h+g^2/b>0$. The projection lower is

$$
fQ_\varepsilon f^*\ge-\theta_\varepsilon^2/b
+(\gamma\zeta+\omega\theta_\varepsilon)^2/a.
$$

Here $\theta_\varepsilon=fQ_\varepsilon y_0^*$ and $\zeta=fPw^*$.
The original theta interval is enlarged by
$\varepsilon\sqrt{A r_{0,\rm up}}$ with a SAVED actual $fPf^*$
UPPER $A$, held fixed for both pricing recipes. Source-energy lower
endpoints never replace this upper. All common identities come from
the pinned reference and rows; interval overlap is only an
inconsistency stop.

After eliminating every real steering coefficient, the comparison
lower is concave in $\theta_\varepsilon$. Its minimum over that
interval is at one of its two endpoints. The saved positive rank gap
$k-g^2/b$ makes the remaining zeta quadratic convex. Outward rational
coefficient bounds and a linear-radius debit reduce its lower to the
exact clipped minimum on the entire saved zeta interval. The clipped
optimizer belongs to a bounding polynomial, not an actual row or
physically legal steering value. Singular $P$ is allowed; failure to
certify the positive denominators or curvature supplies no new floor.

| Quantity | Outward enclosure or directed bound |
|---|---:|
| Transported theta |contained in[−4.331314,−4.325974]|
| Saved zeta |contained in[−8.423007,0.178344]|
| Theta transport error upper |<1.590332e-26|
| Positive-plane energy |>1629.238166139474|
| Fixed comparison curvature |>19789.907603158162|
| Fixed rank gap |>18160.669437018688|
| Previous fixed all-real lower, rounded down |−0.572919446|
| New fixed all-real lower, rounded down |−0.484608778|
| Previous energy-endpoint lower, rounded down |−0.779703304|
| New energy-endpoint lower, rounded down |−0.687084970|

Both floors strictly improve and both exclusion flags remain false.
The four endpoint bounding quadratics select the saved lower zeta
endpoint; this constructs no source realizing that endpoint or a
negative value. A nonpositive lower is inconclusive. Every stronger
previous certificate remains retained. A positive endpoint floor
would exclude only the fixed monotone energy-cap-only pricing family,
with the same row and form, rather than an actual smaller-fee center.

The [exact packet](full-gamma-joint-plane-paid-floor.json.zlib) has
SHA256 `3bdf5f3f81daeeda99c5fa71f2a48ee64c84c09da0f20070b45b0ad781aeec22`;
source SHA256 is
`fb6ae6eb627079e43628481bae7e2b9485077ab22bac7c85b4bcf12994f1b3b9`.
One fixed1024-working/640-export application uses two paid recipes,
two analytic theta endpoints each,12 invoked and66 inherited bindings.
Saved assembly time is<0.133236seconds, excluding serialization.
There is no original producer/fixture replay, new finite/physical
supplier or source/parameter/precision/cutoff scan or retry. This
applies existing indefinite projection and scalar minimization; it
claims no generic theorem, novelty or Lean proof. Original common-
form/domain/source/lift/reference/inverse/weight and actual-symbol
error premises remain. Actual-center/full-border/good-block/cofinal
positivity, Robin and RH remain unresolved.

Reproduction only when explicitly needed, using a fresh output:

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_joint_plane_paid_floor.py --out /tmp/c27-joint-plane-reproduction.json.zlib
```


### A signed source response from saved finite outputs

The archived finite Cauchy block and old approximate jet both contain32768
entries for the SAME original384-bit row0. The old jet A polynomial has a
certified nonzero denominator, so algebra recovers its common approximate
symbol without repeating the source convolution or elementary sine work.
The new direction's complete112features then give one new finite signed
functional. Approximate-symbol and physical-operator errors remain distinct
and are both paid.

Use13 fixed cached Cauchy samples to interpolate
$U(x)=\sum_nz_n/(1+(n/M)x)$, $M=32768$. For rational nodes
$x_i=M/(M+1+\lfloor i(M-1)/12\rfloor)$, the exact Cauchy-kernel
Lagrange remainder supplies the uniform error
$Z\prod_i x_i/(1+x_i)$ on $[0,1/2]$. The old source-power polynomial
is not regenerated. Original real Hurwitz kernels and the saved Abel
phase allowance contract the new coefficients; one cached physical
left-kernel action supplies the mixed correction. Every source interval,
weight, prefactor and separate physical/smooth/forward error is retained.

The fixed1024-working/640-export application gives the following
outward displays; the packet's exact rational bounds and640-bit
interval endpoints remain authoritative.

| Quantity | Enclosure or directed bound |
|---|---:|
| Nonzero old $A_b$ absolute lower |$>1037.446293$|
| Source absolute-sum upper $Z$ |$<2161.865701$|
| Uniform interpolation error $E_U$ |$<0.016074441$|
| Principal-source tail norm error $\delta_v$ |$<0.000017270703$|
| Functional interpolation/symbol error $R_{\rm fun}$ |$<0.000066044862$|
| Oscillatory tail allowance $R_{\rm osc}$ |$<0.000584639478$|
| Physical mixed error allowance |$<4.458569\times10^{-65}$|
| Smooth remainder allowance |$<5.240841\times10^{-12}$|
| Proxy forward allowance |$<3.723887\times10^{-67}$|
| Previous $\zeta$ |contained in$[-8.423006305,0.178343401]$|
| New $\zeta$ |contained in$[-7.591449537,-7.590148168]$|
| Previous fixed all-real lower, rounded down |$-0.484608778$|
| Retained fixed all-real lower, rounded down |$-0.472489776$|
| Previous energy-endpoint all-real lower, rounded down |$-0.687084970$|
| Retained energy-endpoint all-real lower, rounded down |$-0.674685242$|

The actual same-source response interval is strictly narrower, and both
retained paid lower bounds strictly improve. Both exclusion flags remain
false. All four saved bounding quadratics select the new lower zeta
endpoint; this constructs no source realizing that endpoint and supplies
no negative actual paid value. Even after the response is localized to
this narrow negative interval, these particular paid estimates do not
establish positivity. Further sharpening of this scalar alone has not
been shown sufficient.

The application binds16 invoked and66 inherited paths. Stored assembly
time is<1.118949seconds, excluding serialization. The old producers and
accepted fixtures are reused without replay.

The packet retains the recovered approximate symbols, interpolation
coefficients/nodes and physical left-kernel action for later reuse. These
are enclosures of the original source, not arbitrary realizable interval
endpoints. The source and packet SHA256 are respectively
`da93886e4aa00ef58ab4c446486c671996f6f22eaed8d251b7fa8d6d064d94d2` and `2f6fa2f46f91dfa0449d5bec2523cdbd91a0e7587284f8e417ccb5bd2068ee03`.

The one fixed1024/640 application performs one cached-symbol recovery,
one13-node interpolation, one new finite functional, one cached tail
application, one physical left-kernel action and two paid-plane uses.
There is no source/sine/phase/moment/inverse/solve supplier, old producer
or fixture replay, or source/parameter/precision/cutoff scan.
A nonpositive paid lower remains inconclusive. A positive lower excludes
only its specified paid comparison family. Actual smaller-charge center,
full-border/good-block/cofinal positivity, Robin and RH remain unresolved.
This applies existing interpolation/algebra/error estimates; it claims no
new generic theorem, novelty or Lean proof.

[Source](full_gamma_cached_source_response.py) ·
[Exact packet](full-gamma-cached-source-response.json.zlib)


### Bound further refinement of the saved response interval

Keep both saved pricing recipes, the finite parameters, theta endpoints
and directed quadratic coefficients fixed. Allow only a nonempty closed
response subinterval $J\subseteq I$, including a singleton. The issued
lower is the minimum of
$q_{e,J}(z)=A_ez^2+2B_ez+C_e-2R_e\sup_J|z|$ over $z\in J$ and the
two theta endpoints, with $A_e>0,R_e\ge0$.
Its use as an actual comparison lower also requires the true response
$\zeta\in J$. The algebraic algorithm ceiling holds for every such
subinterval, regardless of whether this extra membership is known.

Dropping the nonnegative debit upper-bounds that LOWER certificate.
Convexity bounds each optimistic polynomial on $I$ by
$U_e=\max\{A_ez_-^2+2B_ez_-+C_e,A_ez_+^2+2B_ez_++C_e\}$.
Therefore every such scalar-only refinement, retaining an old lower
$L_0$, has issued certificate at most
$\max\{L_0,\min_eU_e\}$. A negative ceiling excludes positive
certification only by this unchanged algorithm, even with exact
response knowledge. It is not an actual-energy upper or negative row.

The exact-rational application gives the following directed displays.
Lower certificates round down; certificate CEILINGS round up. The
saved fractions remain authoritative.

| Quantity | Fixed recipe | Energy-endpoint recipe |
|---|---:|---:|
| Current retained lower, rounded down |$-0.472489776$|$-0.674685242$|
| Optimistic ceiling at lower theta endpoint |$<-0.472470198$|$<-0.674665225$|
| Optimistic ceiling at upper theta endpoint |$<-0.460850945$|$<-0.663035940$|
| Retained scalar-refinement certificate ceiling |$<-0.472470198$|$<-0.674665225$|

Both recipe flags are true: with the specified coefficients and all
other inputs fixed, no response-only interval refinement can make
these lower certificates positive. This includes a hypothetical exact
response and even waives the nonnegative linear-radius debit. The
result does not exclude actual paid positivity or yield a negative
physical source. It gives no restriction on a new bound that changes
the projection, directed coefficients, prices or joint source
constraints.

The single exact-rational application stores assembly time<0.041523
seconds, excluding serialization. Its eight endpoint evaluations and
two retained caps are new; the original plane and source suppliers
are not repeated.

The exact-rational consumer makes eight new endpoint evaluations and
two recipe caps. It calls no previous plane-floor, physical/source
producer or accepted scientific fixture, and performs no parameter,
source, precision or cutoff scan/retry. Source SHA256 `bc523a5b5392e37f6cbcbb8da27151a39cd3f2e216633e01bb51ff03cd1e1baf`;
packet SHA256 `f4e444527b6e471b710734e702d9d31067a73eac59cab985a358900213729f1a`.

Changing joint source constraints, directed coefficients, theta,
projection or pricing is outside this result. A quantitative positive
lower for the actual source component orthogonal in the same finite
form to the retained timelike/positive plane remains a missing supplier;
no such bound is inferred from a negative certificate ceiling.
This is a classical convexity application, with no generic theorem,
novelty or Lean claim. Actual paid/center/full-border/good-block/cofinal
positivity, Robin and RH remain unresolved.

[Source](full_gamma_response_ceiling.py) ·
[Exact rational packet](full-gamma-response-ceiling.json)


## Test the raw-coordinate cap on the same-source complement

Keep the preceding real form $Q=\gamma P-e_0e_0^*/\delta$, plane
$\operatorname{span}\{y_0,v\}$ and transported response
$\theta=\theta_\varepsilon=Q(f,y_0)$. The accepted factorization
$P=TRT^*$ has an identity tail in $T$, so
$n=(1,\Sigma_{0,\mathrm{tail}}/\Sigma_{00})$ lies in $\ker P$ and
$n[0]=1$. The reference row $\xi=-\delta n$ satisfies
$Q\xi^*=e_0$ and $Q(\xi,\xi)=-\delta$.

Write $\lambda=y_0[0]$ in RAW finite-row units and retain
$\omega=c+g/b$, $a=t+b\omega^2>0$. Its projection off the same plane has

$$
T_\xi=-\delta+\frac{\lambda^2t}{ba},\qquad
\ell_\xi=Q(f,\xi_\perp)
=f[0]+d,\qquad
d=\lambda\left(\frac{t\theta}{ba}-\frac{\gamma\omega\zeta}{a}\right).
$$

The original cap $C\ge|f[0]|$, saved as
`relaxed_negative_coordinate_cap`, uses the SAME raw units as
$\lambda=\alpha_0\rho[0]+\beta_J[0]$; neither is divided by
$\sqrt\delta$. The reference $\xi$ is a row in the full controlled
finite real form domain, not a newly asserted physical source row.
All inherited source/form/physical-error premises remain.

If $I_d$ encloses $d$, then $I_\ell=[-C,C]+I_d$ encloses the actual
$\ell_\xi=Q(f_\perp,\xi_\perp)=f_\perp[0]$. With a certified
$0<T_\xi\le T_{\xi,\mathrm{up}}$, reuse PSD Cauchy to obtain

$$
R\ge r_\xi=\frac{\operatorname{dist}(0,I_\ell)^2}{T_{\xi,\mathrm{up}}}.
$$

The consumer requires a strictly positive directed norm lower before
division. It augments the saved lower for the SAME $F$ by $r_\xi$, then
takes a maximum with the previously retained lower; the residual is
not added to an unrelated earlier bound.

The single1024-working/640-export application gives these outward
ranges:

| Quantity | Enclosure or exact result |
|---|---|
| $T_\xi$ | $[3.460702373153308560,3.460702373153308561]\cdot10^{-31}$ |
| $I_d$ | contained in $[-2.135953295,-2.133377122]\cdot10^{-16}$ |
| $I_\ell$ | contained in $[-6.503116871,2.233786454]\cdot10^{-16}$; exact saved interval contains0 |
| $r_\xi$ | $0$ |
| Fixed retained lower | unchanged, $-0.472489775052793370265\ldots$ |
| Energy-endpoint retained lower | unchanged, $-0.674685241069653953208\ldots$ |

The norm is certified positive; the coordinate cross is not separated
from zero. This estimate gives no positive residual gain. It neither
shows the actual residual vanishes nor supplies an actual energy
upper or a negative physical row.

There is a sharper explanation for this particular interface. The
entire saved $I_d$ lies inside $[-C,C]$, with the accepted
$C>4.367163575891468984\cdot10^{-16}$. For each common realization,
the plane row

$$
f_H=-\frac{\theta}{b}y_0+\frac{\gamma\zeta+\omega\theta}{a}v
$$

has the same two responses, $Q(f_H,y_0)=\theta$ and $P(f_H,w)=\zeta$,
and satisfies $f_H[0]=-d$, $R(f_H)=0$. Thus the two responses and this
absolute coordinate cap ALONE still allow zero residual. This does not
check every earlier scalar cap, a paid sign or generation by the
original source/coupling. A same-source noncancellation relation for
$f[0]+d$, or another residual direction, remains missing.

The new consumer reads four accepted scalar packets, makes one new
coordinate-complement application and two lower additions, and stores
assembly time<0.580244seconds before serialization. It calls no old
plane-floor, source/matrix/inverse/solve/tail supplier or accepted
scientific fixture and performs no parameter/source/precision/cutoff
scan or retry. Source SHA256
`47f34f89cfddfecb581fe626415d303df64de0b5c5004935ebfb0c42285c0a19`;
packet SHA256
`7f6c5d4dc3175e4cd508aadef54d65588bf4c003fddfd01f1c22f8b1748198bb`.
This is an application of existing projection/Cauchy, with no generic
novelty or Lean claim. Actual paid/center/full-border/good-block/cofinal
positivity, Robin and RH remain unresolved.

Reproduction uses the existing pinned Python3.13.12/python-flint0.9.0
environment and an unused output path:

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_coordinate_complement.py --out build/c27-coordinate-complement-reproduction.json.zlib
```

[Source](full_gamma_coordinate_complement.py) ·
[Exact rational/directed interval packet](full-gamma-coordinate-complement.json.zlib)


## Supply a signed raw coordinate from the same cached source

The missing observation can be supplied without regenerating the source.
Keep the original physical map $f=\kappa u_F^*W_FB$ with $B=R_FV$ and
$\kappa=33/32$. Let $e_0$ be the first coordinate of the same18-column
lift, $b_0=Be_0$, and $h_0=(P_0D_{\rm phys}V)e_0$. Then

$$
f[0]=\kappa\langle u_F,b_0\rangle_{W_F}.
$$

This is the raw coordinate used above, with no division by $\sqrt\delta$.
The saved full physical Gram bound supplies
$N_0=\sqrt{G_{00,\mathrm{up}}}\ge\|b_0\|_{W_F}$; its original weight
already belongs to $W_F$ and it contains no $\kappa$ factor. The lift
cap is $L_0=\sqrt{(V^*V)_{00}}$. The accepted forward and mixed-error
certificates are uniform coefficients times a lift norm cap. Their
transport by $L_0/L_b$ therefore applies to this new direction;
it is not an empirical scaling of an error measured on one vector.

Use the stored approximate symbols, source interpolation polynomial
and left physical-kernel action on the first feature column. The
finite principal pairing and analytic tail are evaluated only for
this new column. The mixed term is one new112-component dot product
of the stored left-kernel row with $h_0$, multiplied once by $\kappa$.
No kernel positivity is assumed.

Write $N_{\mathrm{pr}}$ for the saved principal-source norm cap, $\Delta$ for the
saved source-tail interpolation error, $\epsilon_0$ for the transported
physical forward error and $\epsilon_{\rm sm}$ for the saved smooth
source remainder. If $A_{0,k}$ are the new first-column
jet coefficients and $\sigma$ is the original uniform symbol error,
the new symbol-jet allowance is

$$
\beta_0=\operatorname{UP}\frac{2\sqrt w\,\sigma}{\pi\sqrt M}
\sum_{k=0}^{55}|A_{0,k}|\left(\frac M{M+1}\right)^{2k},
\qquad M=32768.
$$

The actual-symbol jet norm is at most $N_0+\epsilon_0$. Consequently
the source/symbol allowance is

$$
\kappa\bigl((N_0+\epsilon_0)\Delta+(N_{\mathrm{pr}}+\Delta)\beta_0\bigr).
$$

Pay the new analytic-tail Abel allowance and transported two-sided
physical mixed error separately, together with
$\kappa\epsilon_{\rm sm}N_0$ and $\kappa N_{\mathrm{pr}}\epsilon_0$.
The approximate symbols are not substituted for actual symbols
without this allowance. The smooth polynomial uses its640-bit export
grid; the uniform symbol error uses its832-bit analytic grid.

Intersect the resulting signed interval $I_0$ with the original
$[-C,C]$. The same-source residual observation is now
$I_{\ell,\mathrm{new}}=(I_0+I_d)\cap I_\ell$, using the SAVED plane
offset and complement norm. An empty intersection stops the application.
With the inherited positive norm certificate, reuse

$$
R\ge r_{\xi,\mathrm{new}}=
\frac{\operatorname{dist}(0,I_{\ell,\mathrm{new}})^2}
{T_{\xi,\mathrm{up}}}.
$$

As before, add this only to the lower of the SAME $F$, then take the
maximum with the previously retained lower. Signed-coordinate
separation, residual gain and paid positivity are separate tests.


The single1024-working/640-export application gives these outward
results; exact saved fractions remain authoritative:

| Quantity | Directed range or lower certificate |
|---|---|
| Raw $f[0]$ | $[2.200218370,2.201145838]\cdot10^{-16}$ |
| New $f_\perp[0]$ enclosure | $[6.426507665,6.776871553]\cdot10^{-18}$ |
| Residual lower $r_{\xi,\mathrm{new}}$ | $\ge0.0001193399383345589672>0$ |
| Fixed retained lower | $\ge-0.472370436$; exact saved lower is still negative |
| Energy-endpoint retained lower | $\ge-0.674565902$; exact saved lower is still negative |

Both paid recipes improve by the same saved residual lower, but neither
certificate is positive. Under the inherited common-source/form and
physical-error premises, this additional signed observation rules out
zero residual for the actual row in this direction. It does not assert
that every row satisfying the older observations has positive residual;
the older zero-residual plane construction remains valid for those
older observations alone. A negative retained lower does not establish
negative actual paid energy.

The source and packet retain SHA256
`a8f9133a4dd8a574fc9a3b7db27cfc1e2d9722764b444fc31cdd636fcd653a12`
and
`12cd1ff9de20ab9083d27ddae812d415ec37c95520de619a6f9813304c2081a5`.
Assembly time is<1.073261seconds, excluding serialization. The consumer
reads ten accepted packets, makes one new finite principal functional,
one new cached analytic-tail application and one new short mixed dot.
It makes no source recovery/interpolation, physical left-kernel action,
quadratic/matrix/inverse/solve/moment/sine/phase supplier or plane-floor
call. No old producer or accepted scientific fixture is replayed;
there is no parameter/source/precision/cutoff scan or retry.

This is a new source-specific application of the existing projection,
error transport and PSD Cauchy argument, with no new generic theorem,
literature novelty or Lean claim. Actual paid/center/full-border/good-block
and cofinal positivity, Robin and RH remain unresolved.


Reproduction uses the existing pinned Python3.13.12/python-flint0.9.0
environment and an unused output path:

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_signed_coordinate.py --out build/c27-signed-coordinate-reproduction.json.zlib
```

[Source](full_gamma_signed_coordinate.py) ·
[Directed interval/rational packet](full-gamma-signed-coordinate.json.zlib)


## Replace one coordinate by the complete same-source far row

Further precision on a single coordinate cannot supply arbitrary
residual gain. Keep the saved cross interval $I_\chi$ and positive
norm interval $[T_{\xi,\mathrm{lo}},T_{\xi,\mathrm{up}}]$. Set
$U_\chi=\max_{x\in I_\chi}|x|$. For any justified refinement $J$ and
norm upper $T_{\mathrm{up,new}}$ enclosing the SAME actual quantities,

$$
\frac{\operatorname{dist}(0,J)^2}{T_{\mathrm{up,new}}}
\le C_\xi:=\frac{U_\chi^2}{T_{\xi,\mathrm{lo}}}.
$$

Indeed $\operatorname{dist}(0,J)\le|f_\perp[0]|\le U_\chi$ and
$T_{\mathrm{up,new}}\ge T_\xi\ge T_{\xi,\mathrm{lo}}>0$.
With the SAME saved plane lower $L_F$ and all other inputs fixed,
the retained LOWER certificate can therefore be no greater than
$\max\{L_{\mathrm{previous}},L_F+C_\xi\}$. This bounds the strength
of this one-coordinate certificate; it is not an upper bound on the
actual residual or paid energy. A negative ceiling rules out positive
certification by that refinement alone, including exact knowledge of
the coordinate and its dual norm.

The full fixed row is a different observation. Keep
$f=\kappa u_F^*W_F B$ in its original18-coordinate lift. The first
coordinate is already saved. Each remaining original basis column is
now supplied once from the SAME accepted source caches, yielding a
joint containing box for $f$. Completeness here means all18entries in
this fixed reduced lift; the unbounded support/source family remains
separate. This constructs no new source and makes
no assertion that independently chosen box endpoints are jointly
attainable. Shared source and error terms do not invalidate a containing
box; evaluating the finite expressions on that box remains conservative.

For each $j=1,\ldots,17$, use the saved feature column $h_j$, lift norm
cap $L_j=\sqrt{(V^*V)_{jj}}$, and physical norm cap

$$
N_j=\sqrt{(G_{\mathrm{far,up}})_{jj}}\ge\|Be_j\|_{W_F}.
$$

The accepted FULL physical Gram Loewner majorant supplies each diagonal
bound. Its original weight is already included, with no factor of
$\kappa$. Reuse the established symbol/source/forward/smooth/mixed
allowances on the new columns, including each new Abel allowance.
The uniform physical errors scale by the SAME $L_j/L_b$ contract.
The stored left-kernel row supplies all17 mixed centers in one matrix
product. The masked kernel is not assumed PSD. Coordinate0 is reused
without another finite or analytic-tail evaluation.

Let $R_{\mathrm{ref}}$ be the saved positive reference inverse and use
the original elimination map $T$, whose first row is
$-\Sigma_{0,\mathrm{tail}}/\Sigma_{00}$ and whose tail is the identity.
Keep

$$
P=TR_{\mathrm{ref}}T^*,\quad
\delta=-\Sigma_{00}>0,\quad
Q=Q_\varepsilon=\gamma P-e_0e_0^*/\delta\le Q_{\mathrm{true}}.
$$

This is the SAME lower comparison form, rather than a newly computed
inverse. The finite rows $\rho,\beta_J$, reference $\alpha_0$ and
$y_0=\alpha_0\rho+\beta_J$ remain unchanged. With $z=fT$, one new
reference action $R_{\mathrm{ref}}z^*$ supplies three contractions:

$$
\begin{aligned}
p_{ff}&=zR_{\mathrm{ref}}z^*,\\
p_{fy}&=(y_0T)R_{\mathrm{ref}}z^*,\\
p_{f\rho}&=(\rho T)R_{\mathrm{ref}}z^*.
\end{aligned}
$$

They give the joint same-source observations

$$
\begin{aligned}
q_{ff}&=Q(f,f)=\gamma p_{ff}-f[0]^2/\delta,\\
\theta_*&=Q(f,y_0)=\gamma p_{fy}-f[0]y_0[0]/\delta,\\
r_*&=Q(f,\rho)=\gamma p_{f\rho}-f[0]\rho[0]/\delta.
\end{aligned}
$$

These use the raw coordinate with its original $\kappa$ factor.
The new $\theta_*$ encloses the SAME actual transported response;
it is supplied directly from the full row and comparison form.
Saved-interval overlap checks detect inconsistencies. The shared
source/form construction supplies the identities; overlap alone does
not establish them.

Reuse $b=-Q(y_0,y_0)>0$, $g=Q(\rho,y_0)$, $h=Q(\rho,\rho)$ and
$a=h+g^2/b>0$. With $s_*=r_*+(g/b)\theta_*$, the complete positive
orthogonal residual is

$$
R_{\mathrm{full}}=q_{ff}+\theta_*^2/b-s_*^2/a\ge0.
$$

Its directed lower is combined with the saved coordinate residual
lower by a MAXIMUM, since both bound the SAME residual. Only afterward
is it added to the cached SAME-$F$ plane lower. The two residual
lower bounds are never added to each other.

There is also a direct route retaining the joint row throughout the
all-real steering completion. For $u=\alpha-\alpha_0$, the SAME
comparison quadratic is

$$
\begin{aligned}
k\alpha^2+D+Q(y_0+f+u\rho,y_0+f+u\rho)
={}&(k+h)u^2+2A_*u+C_*,\\
A_*={}&k\alpha_0+g+r_*,\\
C_*={}&k\alpha_0^2+D-b+2\theta_*+q_{ff}.
\end{aligned}
$$

Under $k+h>0$, its minimum over all real steering coefficients is
$C_*-A_*^2/(k+h)$. A directed interval lower gives another valid
lower for the original paid comparison. Keep the maximum of this
new direct lower, the complete-residual plane lower and the previous
retained lower. They are alternative lower estimates, not summands.
An upper endpoint for the minimum of this LOWER comparison form is
not an upper bound for the actual paid energy.

The single1024-working/640-export application separates three outcomes:
the saved-coordinate ceiling is negative for both retained recipes,
the complete fixed row is supplied, and this rectangular row enclosure
does not improve either paid lower. Exact saved fractions are
authoritative; the following displayed intervals are outward.

| Quantity | Directed range or lower certificate |
|---|---|
| Coordinate-only residual-certificate ceiling $C_\xi$ | $[0.000132707130189352722726116,\ 0.000132707130189352722726117]$ |
| Fixed coordinate-only retained-certificate ceiling | $[-0.472357067922604017542290,\ -0.472357067922604017542289]$ |
| Energy-endpoint coordinate-only retained-certificate ceiling | $[-0.674552533939464600486126,\ -0.674552533939464600486125]$ |
| Complete-row residual interval | contained in $[-7.352650,\ 7.352650]\cdot10^{194}$ |
| Retained residual lower | unchanged, $\ge0.0001193399383345589672$ |
| Fixed direct all-real lower | exact saved lower lies in $[-1.393700,\ -1.393699]\cdot10^{193}$ |
| Energy-endpoint direct all-real lower | exact saved lower lies in $[-1.398417,\ -1.398416]\cdot10^{193}$ |
| Fixed retained paid lower | unchanged; exact lower lies in $[-0.472370435114458811297816,\ -0.472370435114458811297815]$ |
| Energy-endpoint retained paid lower | unchanged; exact lower lies in $[-0.674565901131319394241652,\ -0.674565901131319394241651]$ |

Thus exact refinement of this one coordinate and its dual norm cannot
certify a positive paid lower while the other saved inputs are fixed.
This conclusion concerns the specified lower algorithm, not the actual
residual or paid energy. It does not exclude another joint observation
or a stronger plane lower.

The rectangular joint row remains a valid containing set, but it loses
useful accuracy under the reference contractions. All17 new raw
coordinate widths are at most $1.281033\cdot10^{-19}$. Their stored
finite-evaluation widths are below $2.490399\cdot10^{-43}$, whereas the
separately paid source and Abel allowances are on the $10^{-22}$ to
$10^{-20}$ scales. The resulting positive-reference quadratic interval
is contained in $[-5.871461,\ 5.871461]\cdot10^{96}$; even its sign is
not resolved by interval evaluation. The exact quadratic is nonnegative
under the inherited reference-form premises. Its wide interval does not
contradict that fact.

This diagnoses failure of the present rectangular error transport, not
negative actual energy or failure of the common source. Merely reducing
arithmetic roundoff while keeping these analytic error allowances and
the same rectangular contraction does not remove those allowances.
A stronger use of the row would need a justified estimate preserving
shared-source relations through the reference metric, or another
independently stronger bound. No such improvement is certified here,
and no retry or direction/precision/parameter scan is performed.

The source and packet retain SHA256
`74e900b9d40834c07bc2dc4fd7b37e1bb025542933a1c48073d08a8265bffc41`
and
`00b9e4bdd8a99da4d6fdafc696d313bd43e4060e1ba7d8cec316f7cc5bcb7707`.
The packet is42121bytes; assembly time is<9.746238seconds, excluding
serialization. It reads16 immediate accepted packets, reuses coordinate0,
and makes17 new finite functionals,17 new cached analytic-tail
applications, one17-output mixed product, one reference action, three
reference contractions and two direct minima. It makes no source
recovery/interpolation, physical left-kernel action, old reference
quadratic, inverse, solve or old plane-floor call. Accepted producers,
scientific helpers' tests and fixtures are not replayed.

This is a source-specific application of existing interval enclosure,
projection and scalar square completion, with no generic theorem,
literature novelty or Lean proof claim. All inherited common-source,
physical-error, reference-form and paid-comparison premises remain.
Actual paid/center/full-border/good-block and cofinal positivity, Robin
and RH remain unresolved.


Reproduction uses the existing pinned Python3.13.12/python-flint0.9.0 environment and an unused output path:

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_joint_far_row.py --out build/c27-joint-far-row-reproduction.json.zlib
```

[Source](full_gamma_joint_far_row.py) · [Directed interval/rational packet](full-gamma-joint-far-row.json.zlib)


#### Bound the complete residual before another source query

The full-row rectangular enclosure does not require another source query
to obtain a sharper residual UPPER. Reuse the SAME real comparison form,
raw coordinate, accurate transported response and positive-plane response:

$$
Q=\gamma P-e_0e_0^*/\delta,\quad
R=\gamma P(f,f)-f[0]^2/\delta+\theta^2/b-s^2/a,\quad
s=\gamma\zeta+\omega\theta,\quad \omega=\chi+g/b.
$$

Here $\chi$ denotes the saved coordinate ratio, $b,a,\delta>0$, and the saved scalar
$A_{\rm up}$ bounds precisely $P(f,f)$ for the same fixed18-coordinate
row. The original raw $\kappa=33/32$ normalization is already in $f$;
no additional $\kappa$ or pivot normalization is inserted.
Let $I_0,I_\theta,I_s$ enclose $f[0],\theta,s$ for one common realization.
Set $d_0=\operatorname{dist}(0,I_0)$,
$d_s=\operatorname{dist}(0,I_s)$ and
$m_\theta=\max_{x\in I_\theta}|x|$. Directed denominator bounds give

$$
0\le R\le U_R:=
\gamma A_{\rm up}
-\frac{d_0^2}{\delta_{\rm up}}
+\frac{m_\theta^2}{b_{\rm lo}}
-\frac{d_s^2}{a_{\rm up}}.
$$

Every termwise inequality holds for the actual common realization.
Correlations need not be discarded from the model or asserted to attain
separate interval endpoints. A negative $U_R$, or an upper below the
accepted positive residual lower, is inconsistent and stops the
application. The smaller of $U_R$ and the saved full-row upper remains
an upper for that SAME residual.

This upper can screen a different question from actual paid-energy
sign. Keep a saved lower $L_F$ of the SAME $F$, and an already retained
paid lower $L_{\rm previous}$, fixed. Any valid future lower $r$ of
this residual satisfies $r\le R\le U_R$. Therefore

$$
\max\{L_{\rm previous},L_F+r\}
\le C_R:=\max\{L_{\rm previous},L_F+U_R\}.
$$

A negative $C_R$ prevents positive certification by refining this
residual alone, even with exact knowledge of $R$. A nonnegative $C_R$
is inconclusive. The statement does not bound actual paid energy from
above and does not exclude a stronger joint $F$ estimate or changed
recipe. In particular, $L_F$ is a LOWER for $F$, so substituting it
into an actual-energy upper would reverse the role of the available
information.


The saved-scalar application gives the outward enclosure

$$
0.0001193399383345589672\le R\le
2.364814143855336105831392.
$$

The upper replaces the previous rectangular upper of approximately
$7.35265\cdot10^{194}$. This is an upper-enclosure improvement, not
an improvement of a residual or paid LOWER certificate. With the SAME
saved $L_F$ and previous paid lower, the resulting fixed and
energy-endpoint certificate ceilings respectively lie in

$$
\begin{aligned}
C_{R,\rm fixed}&\in
[1.892324368802542735566375,\ 1.892324368802542735566376],\\
C_{R,\rm endpoint}&\in
[1.690128902785682152622540,\ 1.690128902785682152622541].
\end{aligned}
$$

Both are positive, so this diagnostic does not rule out further
residual estimation. It also certifies no positive paid lower. It uses
only saved scalars; no physical source, feature, matrix, reference,
inverse, solve, analytic tail or previous producer is reevaluated.

#### One common-source covector before bounding errors

Fix exactly one trial by a prescribed rule. Let $e$ be the exact
dyadic midpoint of the saved complete18-coordinate row. Let
$\Sigma_0$ be the exact symmetric midpoint of the original saved core,
and construct $T_0,\delta_0$ from it. Keep the same exact positive
$R_{\rm ref}$ and exact $\gamma$. Define

$$
Q_0=\gamma T_0R_{\rm ref}T_0^*-e_0e_0^*/\delta_0,
\qquad c=Q_0e^*,\qquad z=eT_0.
$$

One exact reference action constructs $c$. This uses $Q_0$ only to
choose a fixed exact covector. It neither replaces the actual $Q$ by
$Q_0$ nor asserts $Q_0\preceq Q_{\rm true}$.
There is no inverse or solve, alternative direction, parameter scan,
precision scan or cutoff change.

The abstract actual-$Q$ Riesz representer of $x\mapsto xc$ is
$u=(Q^{-1}c)^*$. It is used in the proof, not computed. The original
form has inertia one negative and17positive directions, since

$$
Q(x,x)=\gamma(x_{\rm tail}+x[0]t)R_{\rm ref}
(x_{\rm tail}+x[0]t)^*-x[0]^2/\delta,
$$

where $t$ is the first row of $T$. Removing the timelike $y_0$ and
orthogonal spacelike $v=\rho+(g/b)y_0$ leaves a positive definite
space. Define

$$
\begin{aligned}
d_y&=y_0c,& d_v&=\rho c+(g/b)d_y,\\
L_c&=fc+(\theta/b)d_y-(s/a)d_v,\\
D_c&=c^*Q^{-1}c+d_y^2/b-d_v^2/a\ge0.
\end{aligned}
$$

These are the cross and squared norm of the representer projected
to the SAME positive residual space. The Hilbert variational identity
therefore gives, for every real $\tau$,

$$
R\ge2\tau L_c-\tau^2D_c.
$$

Thus an enclosure $L_c\in I_c$ and positive upper $D_c\le U_D$
supply

$$
R\ge r_c:=\frac{\operatorname{dist}(0,I_c)^2}{U_D}.
$$

No positive numerical LOWER for $D_c$ is required. If $U_D<0$, the
inputs contradict structural nonnegativity and the application stops.
If $U_D=0$, then $D_c=L_c=0$ and there is no gain, so no physical
functional need be queried. An interval $I_c$ containing0 gives
$r_c=0$.

The actual-core dependence in the dual norm has an explicit block
formula. Because $c_{\rm tail}=\gamma R_{\rm ref}z^*$,

$$
c^*Q^{-1}c=\gamma zR_{\rm ref}z^*
-\delta(c[0]-t c_{\rm tail})^2.
$$

This supplies no new inverse. Anchor its cancellation exactly at

$$
q_{\rm inv,0}=ec,\qquad
k_0=c[0]-t_0c_{\rm tail}=-e[0]/\delta_0,\qquad
\Delta k=-(t-t_0)c_{\rm tail}.
$$

Then the interval implementation uses

$$
c^*Q^{-1}c=q_{\rm inv,0}
-(\delta-\delta_0)k_0^2
-\delta\bigl(2k_0\Delta k+(\Delta k)^2\bigr).
$$

Finite row actions are also anchored at exact midpoint products, with
all original interval variations retained. The actual $g/b$ remains
in $d_v$; neither the plane nor its readouts are replaced by midpoint
values. The only source-dependent quantities reused in this profile
are the saved accurate $\theta$ and $s$.

For the one new physical quantity $fc$, use linearity to form $Hc$
from the saved112-by18feature map BEFORE estimating any errors.
The complete physical contract is applied to this one aggregate.
A safe lift cap is

$$
L_c^{\rm lift}=\sum_{j=0}^{17}|c_j|
\sqrt{(V^*V)_{jj,\rm up}}\ge\|Vc\|.
$$

The existing uniform sharp-forward and raw mixed allowances scale
by this SAME lift cap divided by their saved positive reference cap;
the raw mixed allowance already includes $\kappa$. One new masked
physical quadratic is evaluated in the original weighted kernel.
Its broad Gram-error cap and its full directional-jet cap are both
paid; take their minimum. The latter retains the omitted signed
antidiagonal moments, anchor bound, moment and prefactor rounding,
and actual-symbol forward allowance. The full high-support moment
Gram does not depend on whether low physical support ends at1280or2048.
No extra $\kappa$ or weight is inserted in this physical norm.

The existing finite functional, cached analytic tail and saved
physical left-kernel action then evaluate this aggregate once.
Original source, approximate-symbol, forward, smooth, raw-mixed and
new Abel allowances all remain. Common source permits this direct
linear functional; it does not authorize guessing correlations from
the previously saved marginal coordinate boxes or discarding any
error payment. Accepted original-coordinate suppliers and their
fixtures are not rerun.

Retain $\max\{r_{\rm old},r_c\}$, check it does not exceed the
saved $U_R$, and add it only to each saved SAME-$F$ lower before taking
the maximum with the previous paid certificate. The two lower bounds
for the SAME residual are never added. A positive output, if obtained,
would concern only the fixed-source paid comparisons under their
inherited premises. It would not settle actual-center/full-border,
cofinal, Robin or RH positivity.


The variational identity is literature-attested: Kurz, Pauly,
Praetorius, Repin and Sebastian, *Functional a posteriori error estimates
for boundary element methods*, Theorem4 proof, lower-bound paragraph,
[DOI10.1007/s00211-021-01188-6](https://doi.org/10.1007/s00211-021-01188-6).
Only its stated Hilbert-space identity is used; no PDE-specific
regularity, approximation or divergence-free conclusion is imported.
The block elimination, Riesz projection and Cauchy steps are standard.
The result here is their source-specific application with the original
physical-error contracts, without a new generic theorem, originality
or Lean-proof claim.


The single fixed1024-working/640-export application supplies the
following outward enclosures. Saved exact fractions and grid endpoints
are authoritative.

| Quantity | Directed range or lower certificate |
|---|---|
| Projected dual norm $D_c$ | $[0.000473029422163194849593401045,\ 0.000473029422163194849593401047]$ |
| New actual aggregate $fc$ | $[-0.062357225190587662283111460115,\ -0.062321776681492458365870736146]$ |
| Projected cross $L_c$ | $[0.000393412340550746747861519848,\ 0.000552651839094430769431308305]$ |
| New retained residual lower $r_c$ | exact lower in $[0.000327195862341560763228991826,\ 0.000327195862341560763228991827]$ |
| Previous retained residual lower | exact lower in $[0.000119339938334558967200315853,\ 0.000119339938334558967200315854]$ |
| Fixed retained paid lower | exact lower in $[-0.472162579190451809501787039320,\ -0.472162579190451809501787039319]$ |
| Energy-endpoint retained paid lower | exact lower in $[-0.674358045207312392445622899610,\ -0.674358045207312392445622899609]$ |

The residual lower increases by more than
$0.000207855924007001796028675973$; both retained paid lower
certificates improve by the same amount. They remain below0, so
neither recipe is certified positive. This does not establish a
negative actual paid energy.

The aggregate physical-norm-squared upper is at most
$0.003645779633808749843771299565$, using the fully paid directional
jet bound instead of the broad Gram cap, which exceeds
$0.089001135794317123782524109570$. These are physical norm bounds,
not $D_c$ or bounds on the source residual $R$.
The new functional's complete explicit allowance is at most
$0.000017724254513063006656728984$; its new Abel part alone is at most
$0.000016648856539851024265837724$.
The projected cross additionally retains the saved plane-response
uncertainties. The successful aggregation therefore resolves part of
the rectangular transport loss while retaining finite, nonzero
analytic payments.

The new source and packet have SHA256
`2379245eacc259969067e7940fa8c1a9b9f7050b0ca901f204fad3d3af1b8cef`
and
`c0aa04f3fbf6d523d0846003a8d2b8385a3d6d3fc891a882d14953a9e88f27c7`.
The packet is48397bytes and assembly takes<0.921542seconds, excluding
serialization. It makes one exact reference action, one core-variation
profile, one aggregated feature column, one masked physical quadratic,
one directional-norm application, one finite functional, one cached
tail and one mixed dot. Original-coordinate queries, new inverse or
solve, old scientific producers and fixture replays are all0.
There is no alternative direction, source, precision, parameter,
cutoff scan or retry.

The saved-scalar envelope source and packet separately retain SHA256
`c9430f825ffa6cb5a745d22e971e2c16c5a53894ee5a4f894a177c7e31180ece`
and
`7ca0ccd8a23dee5a0a60b793d21c07d287652d017afa387e28eb3a6243122d94`.
That packet is13420bytes and assembly takes<0.251660seconds,
excluding serialization. It supplies one upper and two certificate
ceilings without any physical response or old scientific replay.

These are conditional, fixed-source applications of existing methods.
All inherited source, physical-error, reference-form and paid-comparison
premises remain. Actual paid/center/full-border/good-block and cofinal
positivity, Robin and RH remain unresolved.


Reproduction uses the existing pinned Python3.13.12/python-flint0.9.0 environment and unused output paths:

```sh
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_residual_envelope.py --out build/c27-residual-envelope-reproduction.json.zlib
PYTHONDONTWRITEBYTECODE=1 UV_OFFLINE=1 uv run --no-project --python 3.13.12 --with python-flint==0.9.0 python docs/reports/weil-c27/full_gamma_common_trial_residual.py --out build/c27-common-trial-residual-reproduction.json.zlib
```

[Saved-scalar source](full_gamma_residual_envelope.py) · [Saved-scalar packet](full-gamma-residual-envelope.json.zlib) · [Common-trial source](full_gamma_common_trial_residual.py) · [Common-trial packet](full-gamma-common-trial-residual.json.zlib)
