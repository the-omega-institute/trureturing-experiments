# Hardin A185766 recurrence certificate

`check.py` checks the empirical order-14 recurrence in [OEIS A185766](https://oeis.org/A185766), using exact Python integers and the standard library. The source definition is the number of labeled `(n+1)`-by-7 binary arrays whose horizontally or vertically neighboring contiguous 2-by-2 blocks have different ordinary main-diagonal traces; the parent table is [A185769](https://oeis.org/A185769). Boundaries are open and arrays are not identified under symmetries.

Run `python3 check.py`. To use an offline copy of the source data, run `python3 check.py --b-file /path/to/b185766.txt`. The source table is downloaded from `https://oeis.org/A185766/b185766.txt`; its SHA-256 is printed in `output.txt`. OEIS data is attributed to R. H. Hardin and distributed by the OEIS Foundation under [the OEIS licence agreements](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).

The checker reproduces the first 200 terms and all three displayed arrays. Independent tuple-row and bit-row implementations agree on every state and edge. The bit implementation directly counts all admissible two-row and three-row arrays. The recurrence coefficients are supplied by the source; no recurrence fitting occurs.

Let S be the horizontally admissible row pairs, and put an edge `(r,s) -> (s,t)` when all corresponding traces differ vertically. There are 1,736 states and 10,592 edges. Let T sum a vector over outgoing neighbors, and let 1 be the constant vector. The checker verifies, coordinatewise and by both direct and Horner polynomial evaluation, `T^6 q(T)1 = 0`, together with the scalar sums `1^T T^4 q(T)1 = 1^T T^5 q(T)1 = 0`, where `q(X) = X^14 - sum(c_j X^(14-j))` for the printed coefficients.

The intended mathematical argument is unbounded: a bijection of literal arrays with overlapping-row-pair walks gives `a(n) = 1^T T^(n-1)1`; the scalar identities handle n=19,20, and multiplication of the vector identity by every further power of T handles n>=21. The finite identities are helpers for this argument. The n=18 residual is -7776, showing why a stronger threshold would fail.

These are finite arithmetic checks supporting a preregistration. They are not a Lean proof or a claim of kernel verification. The checker source is covered by the repository's Apache-2.0 licence; source sequence data remains subject to the OEIS terms.
