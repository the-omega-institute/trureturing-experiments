# Duro Question 4.17 exact checks

`check.py` uses only Python integer arithmetic and `fractions.Fraction`. It computes complete character table values by the Murnaghan–Nakayama recurrence for every partition pair through `n = 8`, reproducing the source/DHLO values of `D(n)` (including `D(8) = 1`). It then checks the first four terms of the proposed recurrence from `(r,s) = (12,9)`, the triangular and Löschian identities, the exclusion from `B`, and the exact exterior-power hook value on a double transposition.

The first family member is `n = 82`; its partition number is printed by the checker, so a full direct `D(82)` table is outside this bounded fidelity check. The structural family identities and the hook cancellation are checked exactly for `82`, `2705`, `91810`, and `3118757`; the first member also evaluates the binomial coefficient formula directly, while the larger members use the exact vanishing factor before constructing huge binomial integers.

Run:

```text
python3 check.py
```

The stored output ends with `ALL_CHECKS_OK` and the command exits with status 0.
