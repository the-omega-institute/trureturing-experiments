# Six-column trace avoidance (OEIS A185765)

Run `python3 check.py`. Python 3 and its standard library suffice. The stored
reading is `output.txt`; a successful run exits 0 and ends with `ALL_CHECKS_OK`.

The program counts labeled `(n+1) × 6` binary arrays. The trace of the block
at `(i,j)` is the ordinary diagonal sum `A[i][j] + A[i+1][j+1]`. Traces at
horizontal and vertical neighboring block origins must differ. There are no
periodic boundaries or identifications by symmetry.

The source is [R. H. Hardin, OEIS A185765](https://oeis.org/A185765), revision 9,
submitted February 3, 2011, current revision October 3, 2025. The program embeds
the complete [200-term b-file](https://oeis.org/A185765/b185765.txt), which is
source data rather than recurrence-generated test data. Attribution and reuse
terms are the [OEIS End-User License Agreement](https://oeis.org/LICENSE).

Independent literal enumeration checks all 4096 two-row and 262144 three-row
arrays, and all four printed examples. Tuple rows and an independently coded
bit-row oracle agree on every state and edge of the transfer graph. Forward
count propagation matches all 200 source terms. Backward propagation and
integer Horner evaluation verify `T^5 Q(T) 1 = 0` and the two scalar identities
at powers 3 and 4, where

`Q(X) = X^9 - X^8 - 30X^7 - 24X^6 + 95X^5 + 15X^4 - 64X^3 + 116X^2 - 16X - 16`.

These finite identities support an all-length recurrence argument after the
array-to-walk bijection is proved. The experiment is neither a Lean proof of
that bijection nor a Lean proof of the recurrence. It uses exact integers,
without floating point or recurrence fitting. The printed start index is
checked: the recurrence holds through the tested range `13..200`, while its
residual at `n=12` is `-4160`.
