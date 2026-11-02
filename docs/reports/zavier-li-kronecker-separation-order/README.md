# Kronecker separation order check

This exact-arithmetic check reproduces the three-binary-factor witness from Zavier Li, arXiv:2607.12286v1. It stores all logarithms base 2 as integers and uses `fractions.Fraction` for the dual witness; no floating-point value is used as evidence.

The check verifies the diagonal witness

- `H = diag(2^12, 2^-6, 2^-6, 1, 2^-6, 1, 1, 2^6)`, whose identity residual has range 18 and condition number `2^18 = 262144`;
- `G* = diag(2^3, 2, 2, 2^-1, 2, 2^-1, 2^-1, 2^-3)`, whose relative residual has exact eigenvalues `(512, 1/128, 1/128, 2, 1/128, 2, 2, 512)` and range 16, hence condition number `2^16 = 65536`.

For `u,v` with `v < 2u`, it checks the family `log_2 H(z) = 3u(z₁z₂ + z₁z₃ + z₂z₃) + 3v z₁z₂z₃`, the vanishing one-mode partial traces, the exact signed witness giving the additive-factor lower bound `12u + 4v`, and attainment by the local correction `v(z₁+z₂+z₃)`. Tensoring with identity binary factors is checked for `k = 3,…,6`.

Run:

```text
python3 check.py
```

The final line is `ALL_CHECKS_OK`.
