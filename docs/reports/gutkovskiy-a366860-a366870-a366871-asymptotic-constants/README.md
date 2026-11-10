# Gutkovskiy nested floor recurrence constants

`check.py` uses integer arithmetic, `fractions.Fraction` and SymPy exact algebra. Run with Python 3.12 or newer and SymPy 1.14:

```sh
python3 check.py > output.txt
```

For each divisor 3, 4 and 5 the program generates the recurrence through index 1,000,000, checks backward evaluation and the structural inequalities used by the proposed convergence argument, and reports exact rational ratios and error bounds on the finite tail 10,000..1,000,000. It also tests the structural inequalities for divisors 3..40 through index 5,000. The radical formulas, their positive quadratic-root equations and their periodic continued-fraction equations are checked exactly. Rational finite continued-fraction bounds are checked independently.

Fidelity is checked against the initial terms and every entry of Amiram Eldar's 0..10,000 tables linked from [A366860](https://oeis.org/A366860), [A366870](https://oeis.org/A366870) and [A366871](https://oeis.org/A366871). The stored expected digests refer to the comma-separated integer values followed by one newline, without indices. The script's generation code is independently written. No OEIS program code is used. OEIS sequence data is attributed to Ilya Gutkovskiy and the tables to Amiram Eldar; their original [CC BY-SA 4.0 terms](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement) apply to that data.

The source revisions checked are A366860 revision 24, A366870 revision 11 and A366871 revision 11, each carrying a FORMULA conjecture attributed to the entries begun on October 26, 2023. The source pages were read on November 11, 2026. The original program is Apache-2.0 licensed, as specified in the repository LICENSE.

These results verify finite values and algebraic identities. They do not prove convergence at infinity, uniqueness or existence of the infinite recurrence. Those obligations belong to the subsequent mathematical proof. `output.txt` is the successful program output and ends with `ALL_CHECKS_OK`.
