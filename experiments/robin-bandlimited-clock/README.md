# Truncated Fourier-law Robin certificate

The owning project application is
`Library/Estimation/trureturing2026bandlimitedrobin.md`, BFR1–BFR10,
in the trureturing repository. The law is the classical squared Fourier
transform of `sqrt(2) cos(pi v)` on `[-1/2,1/2]`, truncated before its
logarithmic variable is converted to positive multiplicative scales.

```sh
python3 experiments/robin-bandlimited-clock/certificate.py \
  --output experiments/robin-bandlimited-clock/result.json
```

The exact rational comparisons pay the truncated moments, whole-law
normalization and support, the complete complex response at all heights
above `H = 3 * 10^12`, and the original head, infinite tail, directed
transport and elementary allowances. The log enclosure is reused from
`robin-density-transport-3e47/certificate.py`.

For every original integer in the clock slab
`8 * 10^46 < log(n) <= 17 * 10^47`, the complete signed allowance is less
than `81/100`, whereas the reused unrestricted reserve exceeds `81/100`
by more than `1/1000`. The Fourier cancellations, finite zero verification,
five GAK analytic premises and full transport proofs are mathematical
inputs. This program does not numerically evaluate their integrals,
reverify zeta zeros, certify the complete inference in Lean, or prove RH.
No additional zero-density or seven-eighths half-plane input is used.
