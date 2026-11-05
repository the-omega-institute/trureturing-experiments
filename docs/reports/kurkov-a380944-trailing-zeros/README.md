# OEIS A380944 trailing-zero family check

`check.py` implements the OEIS A380944 binary recurrence with `b(v, 0) = 1` and `b(0, k) = 0` for positive `k`, then evaluates `a(v) = b(v, popcount(v))`. It reproduces the 78 published terms, checks the proposed closed form for every `0 ≤ n ≤ 12` and `0 ≤ m ≤ 20`, and checks the finite inclusion-exclusion identity used in the proof outline on the same 273 parameter pairs. All arithmetic is exact; the inclusion-exclusion check uses `fractions.Fraction`.

Run:

```sh
python3 check.py
```

The final line is `ALL_CHECKS_OK` on success. These checks are evidence for the preregistration and are not a Lean proof.
