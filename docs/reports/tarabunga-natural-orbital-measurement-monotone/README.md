# Tarabunga et al. natural-orbital P.4 check

`check.py` uses only `fractions.Fraction` and exact Python integers. It checks the
source fidelity for the even GHZ-like state, the six-mode Gaussian-circuit
branch probabilities and residual parameters, and the strict Rényi-3 average
increase. The final certificate is the integer comparison

`3193^29 * 313^11 * 400^40 < 3364^29 * 484^11 * 343^40`.

The script exits zero and ends with `ALL_CHECKS_OK`.
