# Schulte A194959 dyadic antidiagonal check

Run `python3 check.py`. The saved output is `check.out`.

The program generates the insertion rows directly at zero-based position
`n // 2`. It reproduces all 94 published initial terms and the unweighted
row-sum formula, compares every entry through row 4096 with the published
square-array formula, and checks the conjectured sum on all those rows.
A separate power-enumeration check verifies the bit-based indicator,
including the offset-zero value and the power 2^0. An orientation control
at n=3 distinguishes the reversed antidiagonal in the conjecture from
the forward row. Every calculation uses exact Python integers.

These finite checks test definition matching and small cases. The proposed
unbounded settlement uses a dyadic split and a geometric sum; no finite
positive enumeration is offered as its proof. No Lean proof is included.

Sources: [OEIS A194959](https://oeis.org/A194959), revision 53, with Werner
Schulte's FORMULA contribution and Clark Kimberling's insertion construction;
[OEIS A209229](https://oeis.org/A209229), revision 116. The sequence data and
quoted formulas retain OEIS attribution and its
[CC BY-SA 4.0 terms](https://oeis.org/LICENSE).
