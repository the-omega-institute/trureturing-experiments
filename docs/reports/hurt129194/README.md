# Hurt's A129194 denominator conjecture

`check.py` uses Python's integer arithmetic, `fractions.Fraction` and `math.gcd`; it uses no floating point and no external dependencies. Run `python3 docs/reports/hurt129194/check.py` from the repository root. Its captured standard output is `check.out`; success requires exit 0 and final line `ALL_CHECKS_OK`.

Source: [OEIS A129194](https://oeis.org/A129194), revision 94, offset 0,3. The source NAME is “a(n) = (n/2)^2*(3 - (-1)^n).” Wesley Ivan Hurt's comment (11 July 2016) is “Conjecture: denominator of 4/n - 2/n^2.” The displayed numerical data are attributed to the OEIS Foundation and contributors under [CC BY-SA 4.0](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement). The check program is independently written; no third-party program source is copied.

The fidelity checks reproduce all 53 displayed terms, indexed 0 through 52, and the source's order-six recurrence through index 10000. The settlement checks compare the reduced denominator to the literal source formula at every index 1 through 10000 and at four further large indices, including 10^30 and 10^30+1. They also check the fraction rewrite and gcd cancellation mechanism with exact integers. The source rational expression is undefined at zero, which the program explicitly excludes.

These are finite checks supporting a proposed general proof. They are neither a proof of the universal statement nor a Lean probe. No Lean proof is part of this report.
