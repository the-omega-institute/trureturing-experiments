# Square–decagonal union: OEIS A272266

`check.py` constructs the increasing positive set union of squares and decagonal
numbers directly from Barker's polygonal formula, with repeated values removed.
It uses Python 3.12 or newer and only standard-library integer arithmetic and
`fractions.Fraction`. Run it from any directory with:

```sh
python3 check.py
```

Source: Colin Barker, [OEIS A272266](https://oeis.org/A272266), revision #9,
2016-04-24 08:32:14 EDT; offset 1,2. The entry's definition starts both polygonal
indices at 1. Its FORMULA field calls both the recurrence and generating function
conjectures. The 53 displayed terms and the 1,000-row
[b-file](https://oeis.org/A272266/b272266.txt) supply independent fidelity data.
The b-file is compared through the SHA-256 of its canonical rows, each written as
`n a(n)` followed by a newline:

```text
02b51581077f19062f21833c0a01faf2eeb5e325a2bde90325a5778009186e67
```

The check reproduces all displayed terms and all b-file rows, tests the complete
interleaving inequalities for decagonal indices 2 through 30,000, and compares
all three proposed residue formulas, the literal recurrence, and the formal
rational generating-function coefficients against the first 30,000 semantic
union terms. Sorting the first 30,000 terms of each component gives a complete
prefix: every omitted term exceeds 30,000 squared, while the included squares
already supply 30,000 elements below or at that bound. Polynomial identities are
checked by exact coefficient multiplication. `check.out` is the captured output;
the final line is `ALL_CHECKS_OK`.

These are finite mathematical checks, not a proof at every index or a Lean
verification. The planned all-index argument uses the strict position of each
positive decagonal number after the corresponding odd square and before the next
even square, with the common initial value 1 handled once.

The program was written independently from the mathematical definition. No
third-party program source is copied. OEIS source and data attribution remains
with Colin Barker and the OEIS Foundation; see the
[OEIS licence agreements](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).
