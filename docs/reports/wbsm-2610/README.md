# W-state Bell-state-measurement optimum

`check.py` uses SymPy exact matrices and rational arithmetic to reproduce the
source's `tau_33 = -z` and `tau_13 = (2/3)x` correlations, their attaining value
`Lambda = 13/9`, the standard-frame index-pair upper check, the state-specific
quadratic identity on rational unit vectors, and white-noise scaling. It uses no
floating-point evidence and exits only after printing `ALL_CHECKS_OK`.

Run:

```sh
python3 check.py > check.out
```
