# Fixed unwindowed transfer of the original theta variance gap

This application targets the original even minimal mixed-jump form
with every prime power retained. The conditional lower bound already
owned by [TC12](target-correction.md#reuse-the-positive-blocks-when-a-requested-gap-test-fails)
is $921/2000=0.4605$. The fixed application below targets $15/32=0.46875$,
an increase of $33/4000$. Its reusable Fourier, form-domain,
whole-projector, theta-tail and actual-zero-count arguments are existing
results; no generic theorem or literature-priority claim is made.

## Immutable inputs and fixed comparison

The [program](verified_height_unwindowed_gap.py) reads three accepted
results at experiment commit
`e5c5d8bf2b157e5c94af05ba724e18f4131b4c88`:

- [derivative norms and symbol upper bound](derivative-bandwidth-result.json);
- [128 complete pointwise row pairs](joint-high-floor-result.json);
- [the full exterior potential floor](sharper-exterior-result.json).

It acquires no theta values, operator actions, Grams, moments, zero
ordinates or source data. All comparisons use Python standard-library
exact rational arithmetic on the saved dyadic endpoints. No accepted
producer or fixture is rerun. The construction fixes
$\varepsilon=1/32$, $N=2^{40}$ and $H=3\cdot10^{12}$ and has no scientific
options or scans. It refuses to overwrite an existing result.

Let $J_i$ bound the complete old $W+m(32)s^2$ on cell $i$, and let
$s_{i,-}^2$ bound the same actual $s^2$ there. Existing symbol estimates
give $m(N/2)>40/3$ and $m(32)<22/3$. The new finite question is therefore

$$
J_i+6s_{i,-}^2>31/64\quad(0\le i<128).
$$

The saved cells partition $[0,3/2]$ with width $3/256$. Reflection and
the saved `exterior_W_floor_lower` cover the whole line. The previous
leakage-subtracted high-form floor is not reused at the new bandwidth.

## Fixed result

All 128 row inequalities and every serialized budget check pass. The
smallest row margin is exactly

$$
\frac{721580297331903231342176262114335387}
 {340282366920938463463374607431768211456}>0.
$$

The [result packet](verified-height-unwindowed-gap-result.json) retains
the exact commutator, projector, leakage, tail and contradiction budgets.
No additional numerical supplier is inferred from rounded displays.

## Analytic bridge and remaining boundary

The main repository's [owning proof](https://github.com/the-omega-institute/trureturing/blob/lane/math/rh-c27-propagated-prime-row/Library/Weil/fukushima2011dirichlet.md#a-fixed-unwindowed-transfer-to-the-original-variance-gap)
establishes the specified noncompact paired zero identity on
$W_Nh=U^{-1}p_NUh$ by the existing scaled-cutoff passage in the actual
minimal form. It keeps the original theta normalization and transports
all prime powers. A fixed smoothed trapezoid gives $L_p\le16/15$, $C_p<10$.
The already proved principal commutator coefficient is $64$; the
conservative existing $128$ definition and its consumers remain intact.

The fixed budget checks combine the original derivative packet with
the elementary theta tail estimates in the owning proof. They yield
$K_p^\sharp<9\cdot2^{18}$, the whole-projector form error $<1/80$,
and the actual shifted paired-zero-tail error $<2^{-19}$ using the
published Platt–Trudgian inclusive verified head. The strict ordinate
tail retains all multiplicities and both off-line paired factors.

For a centered unit vector in the low projector, the same-map form
comparison gives $q(W_Nh)<-3/6400$, which contradicts the paid tail
lower bound. Subject to every inherited actual-model, original-domain,
ground, sound-supplier and published finite-height premise, the
conclusion is

$$
D(h)\ge\frac{15}{32}\operatorname{Var}_\nu(h),
\qquad h\in\mathcal F_{\min,\mathrm{even}}.
$$

The [saved result](verified-height-unwindowed-gap-result.json) records
the complete fixed comparison and exact rational budgets. It checks a
new inequality on accepted data, not the analytic premises themselves.
This is a conditional paper deduction, without Lean certification or
an exhaustive novelty finding. The endpoint $1/2$, a cofinal family,
the full Robin inequality and RH remain unresolved. At a fixed verified
height the frequency-separation requirement still limits this supplier.

```sh
python3 docs/reports/theta-mixed-matrix/verified_height_unwindowed_gap.py
```

The program is project-authored and requires only Python's standard
library. An existing saved result is reused without rerunning it.
