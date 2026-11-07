# Finch Conjecture 23: capture offset at k = 2

Source: Steven Finch, *Three Quadratic Recurrences Modulo 3: What We Can Prove and What We Can Guess*, [arXiv:2610.34300v1](https://arxiv.org/abs/2610.34300v1), Conjecture 23, Definitions 15 and 21–22. The source's sequence uses positive indices, arithmetic modulo 3, and a capture run starting at r ≥ k + 1.

Run with Python 3 (standard library only):

```sh
python3 docs/reports/finch-conjecture-23/check.py
```

`check.out` stores the output. The check reproduces all eleven periods in the source's main-text table for k = 2 through 12, and the source's pre-period μ = 5 and capture point z = 11 at k = 3. It also reproduces the offset 6 for the eight source-listed parameters 3, 6, 7, 10, 11, 14, 15 and 19.

At k = 2 the initial state (0, 1) recurs after four transitions, and the recurring scalar block is (0, 1, 1, 2). Every smaller positive period fails on a recurring phase. Thus the least eventual period is 4 and the least pre-period is 0. Here G = 1. The first allowed zero pair is at 5 and 9, giving z = 9 and z − μ = 9, which differs from the claimed 6.

The program uses exact integer arithmetic. Complete order-k state collisions certify infinite periodic tails; scalar prefixes alone are not treated as periodicity evidence. The separate phase test at k = 2 checks the recurrence on every phase of the proposed block. These are computational checks, with no Lean proof or claim of a kernel-checked settlement. No corrected uniform conjecture or infinite counterexample family is asserted.
