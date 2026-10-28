# Quadratic recurrence modulo 3: zero runs

For recurrence (A) in Steven Finch, *Three Quadratic Recurrences Modulo 3: What We Can Prove and What We Can Guess*, arXiv:2610.34300v1, section 6, Remark 20, the universal attainment reading fails at k = 61. The maximum number of gaps in a Definition 15 run is 10, below the Theorem 16 upper bound 15.

Run with Python 3, using only the standard library and exact integer arithmetic:

```sh
python3 docs/reports/quadratic-recurrence-mod3-zero-runs/check.py
```

The stored output is `output.txt`; its final line is `ALL_CHECKS_OK`.

The fidelity check reproduces Example 18's complete 60-term table and zero indices at k = 4. At k = 61, two independent representations agree: the scalar recurrence and a 61-digit ternary state, with the newest digit least significant. The first repeated state has indices 739 and 9685, giving cycle length 8946. This period was already reported in Finch's *Periodicity in Sequences Mod 3* (2002), page 2; discovery of the period is not claimed here.

The 9685 distinct states form a forward invariant set: every state advances to its successor in the list, and the last advances to state 739. None has a 45-digit suffix whose zeros are exactly the positions divisible by four. Every run with at least eleven gaps would produce such a suffix (in reverse order, which preserves its zero mask). Thus the finite certificate addresses the entire infinite canonical orbit, including windows crossing the cycle boundary. The direct check also verifies the ten-gap run at indices 2780 through 2820, with values (0,1,2,2) repeated ten times and a final zero.

This is a reproducible arithmetic certificate and a proposed refutation route. It is not a Lean proof, and it does not refute the upper bound, Conjecture 14, or Conjecture 23. No infinite family of exceptions is claimed.

Sources: [arXiv v1](https://arxiv.org/abs/2610.34300v1), [TeX archive](https://arxiv.org/e-print/2610.34300v1), and [Finch's 2002 note](https://web.archive.org/web/20150911081920id_/http://www.people.fas.harvard.edu/~sfinch/csolve/seqmod3.pdf). The checker is independently written; no source implementation is copied.
