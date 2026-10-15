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
