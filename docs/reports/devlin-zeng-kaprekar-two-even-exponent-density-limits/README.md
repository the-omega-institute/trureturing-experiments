# Devlin–Zeng even-exponent density limits

`check.py` uses exact integer arithmetic and `fractions.Fraction` to reproduce the source's four-digit Kaprekar trajectory and its odd-factor density formula, then evaluates the difference-pair reduction at bases `80`, `320`, and `1280`. It checks the 76 residue classes proposed as an attracted family for exponents divisible by four and the 8 residue classes proposed as a non-fixed periodic obstruction for exponents congruent to two modulo four. The displayed source line `3322 − 2233 = 889` has an arithmetic typo; the script records the correct value `1089` and checks the subsequent exact trajectory.

Run:

```text
python3 check.py
```

The final line must be `ALL_CHECKS_OK`.
