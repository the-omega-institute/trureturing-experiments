# Hardin–Barker A230393: five-colour successor strips

`check.py` checks the literal n by 2 array count in OEIS A230393 and the proposed five-state transfer, using exact integers and SymPy algebra. Run `python3 check.py` with Python 3.12 or newer and SymPy 1.14.0. The checked output is `check.out`; its final line is `ALL_CHECKS_OK`.

The source is [OEIS A230393](https://oeis.org/A230393), revision 12, with offset 1,5. The entry is by R. H. Hardin (2013); Colin Barker's empirical generating function is attributed to 17 March 2018. `b230393.txt` is Hardin's [table for heights 1 through 210](https://oeis.org/A230393/b230393.txt), obtained from OEIS under its [End-User License Agreement](https://oeis.org/LICENSE). Its SHA-256 is `d47859cfd3244ba8f42e364a215c5ddb354d08d9ed479c0b5bd3df7da6ad474f`.

The program exhausts all arrays with upper-left entry zero at heights 1 through 5, checks a separate literal row-pair dynamic program through height 14, and compares the unpruned pending-test transfer and its reflection quotient with every one of the 210 table values. It checks both boundary rules, every quotient transition, the determinant and the relevant adjugate entry of I-xT, the matrix polynomial identity, a signed change of basis for the existing scalar-recurrence API, and the generating-function coefficients and recurrence on the available finite range.

A state is (delta, left-pending, right-pending): delta is the right-minus-left colour modulo 5, and a pending flag means that neither the horizontal neighbour nor the row above has supplied the cell's successor. The next row must discharge both old pending tests. A final row is accepted only when neither flag is pending. The source has no periodic boundary identification and no quotient of arrays by reflections. Reflection is used only to count continuation states; the initial weight is two.

These are finite exact checks supporting a preregistration, not an all-height proof. The intended proof must establish the prefix-state invariant and quotient counting formula for arbitrary height before using the matrix relation. No Lean proof of the assertion or its negation is included here.
