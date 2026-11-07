# Wang--Yuan positive-scattering reshaping check

`check.py` implements the exact two-component specialization of Wang and Yuan's activation formula. Each mode cost is computed from integer supports and normalized coordinate functionals using `fractions.Fraction`; grouping uses exact Kronecker products. The two families use only nonnegative integer factors with overlapping supports and finite scattering, and every grouped decomposition has at least three modes.

For each `k` in `{1, 2, 3}`, the increase family gives `tau = 1` and `tau_grouped = k + 1`; the decrease family gives `tau = k + 2` and `tau_grouped = 2`. The output is stored in `check.out`, and the checker exits zero with `ALL_CHECKS_OK`.
