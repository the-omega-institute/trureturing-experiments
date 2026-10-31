# Wiseman A083368 run-profile first column

`check.py` tests the July 2025 conjecture in [A083368](https://oeis.org/A083368/internal), using its historical highest-changed-position definition, and the first-column observation in [A385817](https://oeis.org/A385817/internal). It uses the low-order-first one-run convention in [A245563](https://oeis.org/A245563/internal); zero has an empty profile.

Run with Python 3 and no third-party packages:

```sh
python3 check.py > check.out
```

The source-fidelity checks reproduce A245563's printed rows 0 through 21, A083368's first data line, and its example for 27 and 28. The program scans every integer below 2^20, independently enumerates the lazy Fibonacci representations by actual weighted values, and compares the complete first-occurrence profile column in that range. It also tests the block identities used by the proposed general proof. All evidence uses exact integer arithmetic. The saved output is `check.out`.

These are finite checks. They establish neither the unbounded identity nor any Lean theorem. The later shortened A083368 NAME differs from the historical statistic already at n=3; that distinction is explicitly tested.

The check program is independently authored for this experiment. Short numerical examples are attributed to the linked OEIS entries; OEIS material is available under its [CC BY-SA 4.0 terms](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).
