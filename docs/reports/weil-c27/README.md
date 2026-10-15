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
