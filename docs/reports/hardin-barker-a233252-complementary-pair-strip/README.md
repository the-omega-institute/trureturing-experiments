# OEIS A233252 complementary-pair strip check

`check.py` is an exact Python 3 fidelity and settlement check for OEIS A233252 (revision 8, 2018-10-10). It enumerates the literal six-color arrays (through canonical complementary-pair classes and all free signs) on the `n × 4` graph for `n = 1, 2, 3`, reproducing `10, 1184, 166400`, and checks the proposed closed form, two growth rates, and recurrence against all 12 displayed OEIS terms. It uses integer arithmetic and `fractions.Fraction`; no floating point is used as evidence.

Run:

```sh
python3 check.py
```

The final output line is `ALL_CHECKS_OK`.
