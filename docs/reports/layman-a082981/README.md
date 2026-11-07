# OEIS A082981 finite check

`check.py` uses only exact Python integers and `fractions.Fraction`. It runs 250 rounds of Layman's simultaneous ratio-determined insertion rule, retaining the first 3000 entries, and checks that the first 200 entries are unchanged in the last round. It matches all 41 terms in OEIS revision 22, matches the exact floor recurrence generator, checks 192 lag-four identities, checks the exact coefficient cycle `(2, 2, 3, 2)`, and verifies the integer matrix identity `B^2 - 6B + I = 0`.

This is a finite consistency check; it is not a proof certificate. The supplied candidate mentioned 248 rounds; with this direct prefix simulation the 200th retained entry changes between rounds 248 and 249, so 250 rounds are used to obtain a 200-term unchanged prefix. This is a finite-check accounting detail and does not settle the infinite assertion.
