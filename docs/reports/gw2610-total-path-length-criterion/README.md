# Total path length criterion: exact finite checks

This report checks the variance identity and the innovation recurrence from
Gałązka--Wdowicka, *Total Path Length in Power-Weight Recursive Trees:
Martingale Limits and Global Fluctuations*, arXiv:2610.00614v1.

`check.py` uses only `fractions.Fraction` and SymPy exact algebra.  It
enumerates every parent sequence through `n = 8` for the integer power weights
`w_n = n`, and through `n = 3` for the exact `c = 1/2`, `gamma = 3/4` family.
It checks the parent-law mass, `eq:v-rec`, the bound
`v_n^2 <= (n - 2)/2`, `eq:VarT-exact`, and the finite truncation identity
`S_N = C_N`.  The latter family is also checked symbolically at exact fourth
powers, where `n p_n = q/2`, giving an unbounded `n p_n` witness.

Run from this directory with:

```sh
python3 check.py
```

The captured output is in `output.txt`; the final line is `ALL_CHECKS_OK`.
