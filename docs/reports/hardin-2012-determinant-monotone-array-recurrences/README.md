# Hardin determinant-monotone binary arrays

`python3 check.py` runs exact integer checks for OEIS A204609, A204610 and
A204611, the widths 2, 3 and 4 of the table A204616. Python's standard library
is sufficient. A successful run exits zero and ends with `ALL_CHECKS_OK`.
`output.txt` is the corresponding complete output.

`source-data.json` records the OEIS entry revisions and the 210-term b-files
linked by each entry. Source text and sequence data are credited to R. H.
Hardin and Colin Barker; the OEIS End-User License Agreement is at
<https://oeis.org/LICENSE>. The source URLs and revision records are included
with the data. The first A204609 generating-function conjecture is checked;
its separate radical formula is outside the target.

The program constructs the last-row/last-minor-vector states directly from
all ordered binary row pairs. There are 8, 18 and 34 states. It prints the
states, initial multiplicities and all transition columns. Independent
exhaustive enumeration using integer-coded rows checks sizes n=1..8, 1..5
and 1..4, respectively. Transfer counts equal all 210 source terms per entry.

For each source recurrence, let p be its monic characteristic polynomial and
b the initial state vector. The program checks T^3 p(T)b = 0 for widths 2
and 3. At width 4 it checks T^3 p(T)b = 20R, TR = R and sum(R) = 0;
the nonzero residual is invisible to the total count. It also checks the
A204609 vector polynomial certificate (I-xT)q = Db and sum(q) = P.
These finite exact identities supply a route to an all-size proof after the
array/transfer correspondence and the invariant propagation are proved.
They are experimental checks, not Lean proofs or completed settlements.
