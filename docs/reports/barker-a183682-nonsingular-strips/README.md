# Hardin–Barker A183682: nonsingular binary strips

Exact finite checks for the number of labelled (n+1) by 3 binary arrays whose consecutive 2 by 2 subblocks have nonzero determinant over the integers. The proposed all-height settlement uses a row-compatibility graph and induction. This report contains consistency evidence, not a proof of the conjecture.

Run with Python 3 and its standard library:

```sh
python3 docs/reports/barker-a183682-nonsingular-strips/check.py
```

The script verifies all 64 ordered row pairs and their two-class quotient, Q squared = Q + 4I, all 200 published terms by independent eight-state propagation, direct labelled-array enumeration for n=1 through 5, the recurrence for n=3 through 200, and coefficients through degree 200 of 4x(3+5x)/(1-x-4x²), using integers and fractions only. The final line of the stored output is ALL_CHECKS_OK. The source recurrence is deliberately not applied at n=2: the literal source-object count at n=0 is 8, whereas the auxiliary active-vertex walk count at length zero is 5. The ordinary generating series has constant coefficient zero.

Source: R. H. Hardin, [OEIS A183682](https://oeis.org/A183682/internal), definition contributed January 6, 2011; Colin Barker's generating-function conjecture, February 27, 2018. The retrieved entry identifies revision 9, February 27, 2018, and retains the labels “Empirical” and “Conjectures.” The radical formula printed in the entry is outside this report's target.

The fixture b183682.txt is the complete [OEIS table by R. H. Hardin](https://oeis.org/A183682/b183682.txt), indices 1 through 200, retrieved November 3, 2026. OEIS sequence data are used under [CC BY-SA 4.0](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement); attribution and the original table are retained here. The check program is original report code and is released under CC BY-SA 4.0.

- check.py SHA-256: a5f0783f8724e0b3409af9a851949d76bf69e7a0340db0bfe7bab4f0430936a5
- b183682.txt SHA-256: 83d9f74426b679f0ed9ec4d9439ddcf0746f8da19b031b9840c91f7ad9d0a53b
- check.out is the exact standard output of the command above, which exited 0.
