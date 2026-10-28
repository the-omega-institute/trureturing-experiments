# OEIS A210196: exact recurrence and coefficient check

Source: [OEIS A210196](https://oeis.org/A210196), revision 21
(2022-02-22), with Clark Kimberling’s defining polynomial pair and Knud
Werner’s FORMULA conjecture (2022-01-10).

Run `python3 check.py`. Python 3’s standard library suffices. The program
uses integer coefficient arrays and the original simultaneous recurrence;
the expected coefficients are computed independently with `math.comb`.
The five printed source rows are transcribed as fidelity fixtures, with
attribution to Clark Kimberling and OEIS. Source data retains the
[OEIS license terms](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).

The stored output `check.out` reproduces all five source rows, checks 5,550
coefficient equalities for rows 1 through 100 (including five indices past
each row), and verifies the source’s alternating row sums for these rows.
The final line is `ALL_CHECKS_OK`. Any failed comparison raises an exception
and produces a nonzero exit status, including when Python assertions are disabled.

This is an exact finite check. It does not prove the universally quantified
coefficient identity. It contains no Lean proof or kernel verification.
