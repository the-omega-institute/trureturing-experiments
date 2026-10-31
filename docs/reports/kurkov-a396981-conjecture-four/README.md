# Kurkov A396981 Conjecture 4: exact finite checks

Source: Mikhail Kurkov, [OEIS A396981](https://oeis.org/A396981/internal),
approved revision 45, June 19, 2026, 09:41:43 EDT; entry authored June 12, 2026.
The source's Python program is attributed to Peter Luschny. The target is the
three-parameter identity explicitly labelled Conjecture 4 in FORMULA.

Run with Python 3 using only its standard library:

```sh
python3 docs/reports/kurkov-a396981-conjecture-four/check.py
```

`check.py` constructs the array using exact integer arithmetic and triangular
block comparisons. It reproduces all 66 published DATA terms and the 72 entries
of the displayed rectangle, and agrees with an independent evaluation of the
source's integer-square-root recurrence on 729 pairs. It checks the proposed
row-difference law on 8,400 pairs and Conjecture 4 on 4,590 triples, including
zero parameters. `output.txt` is its output. Every comparison passed, with final
line `ALL_CHECKS_OK` and exit code 0.

The unused, unlabelled Stirling FORMULA line in the same source is inconsistent
with the defining recurrence at n = 1, k = 3: it gives 6, whereas the recurrence
gives 3. Neither this formula nor a corrected version is used in the checks.

These finite checks establish source fidelity and the stated finite ranges.
They do not prove the identity for unbounded parameters and are not Lean kernel
verification.
