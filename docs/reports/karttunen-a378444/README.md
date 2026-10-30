# Karttunen's A378444 ceiling conjecture

Source: [OEIS A378444](https://oeis.org/A378444), final FORMULA line,
conjectured by Antti Karttunen. The definitions use A083345's reduced rational
numerator and A174273's inverse Möbius transform of A035263.

Run with Python 3 and its standard library:

```sh
python3 docs/reports/karttunen-a378444/check.py
```

The program downloads two OEIS b-files and verifies their pinned SHA-256 digests.
For an offline run, supply `--bfile-dir DIRECTORY` containing `b378444.txt` and
`b174273.txt` with those same digests. `output.txt` contains the deterministic output.

For every integer from 1 through 65537, exact `Fraction` arithmetic agrees with
an independent arithmetic-derivative/gcd numerator computation. Direct divisor
sums match all 65537 published values for both sequences, the proposed ceiling
identity, and the exponent-box product formula, with zero mismatches. Seven
additional finite exponent boxes include the empty box, a power of two with
exponent 1000, and odd and mixed prime powers. The final line is `ALL_CHECKS_OK`.
These finite checks support a preregistration; they do not prove the identity
for all positive integers.

A174273's ancillary formula involving sigma is inconsistent at 3: its defining
transform and b-file give 2, whereas sigma(3) is 4. The program records this
explicitly and uses the defining transform throughout.

The program is original checking code. The externally hosted numeric inputs are
OEIS b-files linked above; no b-file copies are included. Their source and byte
digests are printed in the output. OEIS data are subject to the
[OEIS license agreements](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).
