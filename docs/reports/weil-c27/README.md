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
