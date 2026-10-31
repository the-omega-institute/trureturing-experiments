# Jones–Kinnersley hypercube directional localization

This report checks the conjecture in Question 6.1 of Jones and Kinnersley,
*The directional localization game on graphs*, arXiv:2609.01745v1.

`check.py` uses only exact integer combinatorics. It performs an exhaustive
belief-state game search for dimensions 1, 2, and 3, reproducing the source's
bound `n ≤ dirloc(Q_n) ≤ n+1` and finding the conjectured value `n` in those
three cases. It also checks the proposed general n-cop upper strategy for all
legal replies in dimensions 1 through 5, and checks the paired-history lower
invariant for every closed-neighborhood and projected-probe-set case through
dimension 4. The check is a finite sanity check; the preregistration's
settlement target is the all-dimensional proof.

Run:

```sh
python3 check.py
```

The stored output is `check.out`; the final line is `ALL_CHECKS_OK`.
