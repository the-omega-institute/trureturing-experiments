# Kurkov's two array identities for A006351

The exact check uses the recurrences in OEIS A379459 and A379460 and the
independent Gutkovskiy recurrence in A006351. It reproduces 28 published
antidiagonal entries of each array and 11 published terms of A006351, then
checks both first-column identities for every row from 0 through 30.
It also verifies the source EGF reversion equation through degree 32 using
truncated rational power series. All evidence uses Python integers and
`fractions.Fraction`; no floating point or external package is required.

Run from the repository root:

```sh
python3 docs/reports/kurkov-a006351-two-array-identities/check.py
```

The output is in `output.txt`; its final line is `ALL_CHECKS_OK`.
These are finite fidelity and indexing checks, not a proof for every index.

Sources: [A006351](https://oeis.org/A006351), revision 195;
[A379459](https://oeis.org/A379459), revision 4; and
[A379460](https://oeis.org/A379460), revision 4. Mathematical definitions and
DATA are transcribed; the source's PARI programs are not used or copied.
OEIS content is governed by the [OEIS End-User License](https://oeis.org/LICENSE).
