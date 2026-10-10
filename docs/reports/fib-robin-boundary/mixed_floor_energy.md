# Mixed floor energy and reciprocal-square diagnostics

The [program](mixed_floor_energy.py) evaluates the classical sawtooth
quadratic expression for the fixed mixed-sieve coefficients in
[the source card](https://github.com/the-omega-institute/trureturing/blob/6653df2ce7c0bae76baee5eed1ea88b53ef768e5/Library/ArithSums/srivastav2025sievevaughan.md).
It reuses the coefficient constructors and recorded high-precision means
from [the existing diagnostic](mixed_sieve_response.py); it does not
enumerate that diagnostic's response grid again. All clocks are
$A=m^8$, not actual selected Robin sources.

For $F(y)=\sum_dh(d)\{y/d\}$ and $c=\tfrac12\sum_dh(d)$, the
computed centered periodic energy is

$$
\mathcal E=\frac1{12}\sum_qJ_2(q)
 \left(\sum_{q\mid d}\frac{h(d)}d\right)^2,
\qquad J_2(q)=q^2\prod_{p\mid q}(1-p^{-2}).
$$

The program also samples
$E(y)=(F(y)-y\max(S,0))_+$ at $y=x/p^2$ for every prime
$p^2<x/U$, with five $x$ values per scale:
$A,5A/4,3A/2,7A/4,2A$. Its recorded cost is
$x^{-1}\sum_p\log p\,E(x/p^2)$, using the upper envelope
$\delta_{N,p,2}\le\epsilon\log p$ rather than constructing actual
CA marginal weights. Rational arithmetic determines each fractional
part before multiplication by floating-point coefficients.

| $m$ | Centered energy | Standard deviation / coefficient absolute sum | Square samples |
|---|---:|---:|---:|
| 2 | 0.142745840 | 0.088156521 | 21 |
| 3 | 0.289501871 | 0.044413083 | 53 |
| 4 | 0.538610625 | 0.026470551 | 106 |
| 6 | 1.441464661 | 0.013275072 | 279 |
| 8 | 2.995533052 | 0.007826401 | 567 |
| 12 | 8.926861013 | 0.003832486 | 1583 |
| 16 | 19.911943003 | 0.002255574 | 3330 |

The program evaluates 5,939 square samples; the
[data](mixed_floor_energy.json) record per-clock cost summaries and
selected direct checks. Their computed positive-part costs are zero. All reused means in
these seven finite scales are positive; their drift is retained in
the sampled residual. These observations neither determine its
asymptotic sign nor construct the selected same-source family.

The Jordan-channel expression agrees with direct pairwise gcd
summation at $m=2,3,4$. The smallest sampled response at each scale
is additionally evaluated by the original pair formula at 80 and
120 Decimal digits; their maximum difference is below $4\cdot10^{-79}$,
and the float comparison error is below $1.6\cdot10^{-15}$. Those
checks test these implementations, not an interval or Lean proof.

The period-energy norm is smaller than the coefficient absolute-sum
envelope on these inputs. Transferring it to a uniform supremum over
the actual reciprocal-prime-power samples is still a missing
estimate. A global-period norm cannot be substituted for that
sampling law. No scaling exponent, uniform bound, mean requirement,
MC estimate, Robin improvement or RH conclusion follows from the
finite readings. Only the square layer and five values of $x$ per
scale are sampled; no bound on the complete higher-power component
is claimed.

Run from the private experiment repository:

```sh
python3 docs/reports/fib-robin-boundary/mixed_floor_energy.py
```

The program uses only the Python standard library and the sibling
coefficient module. Its output records the SHA-256 of the reused
mean dataset.
