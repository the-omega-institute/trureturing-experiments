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
