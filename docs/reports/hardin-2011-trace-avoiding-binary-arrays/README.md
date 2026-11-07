# Trace-avoiding binary arrays: OEIS A185761–A185764

These exact checks concern R. H. Hardin's counts of `(n+1) × k` binary arrays whose contiguous `2 × 2` block traces differ from every horizontal and vertical neighboring block trace, for widths `k = 2,3,4,5`. The recurrence and generating-function conjectures are recorded in the accompanying OEIS entry snapshots; Colin Barker supplied the generating functions in 2018.

Run `python3 check.py` with Python 3. The program uses only the standard library and exact integer and rational arithmetic. Its captured output is `output.txt`, ending in `ALL_CHECKS_OK`.

The program checks all 200 published b-file terms for each member by two transfers: full pairs of columns and the compressed terminal trace/current-column state modulo bit complement. The latter has 3, 8, 20 and 52 states, respectively; the unquotiented compressed state counts are 6, 16, 40 and 104. Every state is reachable already from a two-column initial array. Representative independence and transition closure are checked exhaustively. Width 6 has 132 quotient states and is excluded by the 64-state selection bound.

For each conjectured polynomial `q`, the output supplies an integer vector `w = T^b 1` with `q(T)w = 0`, plus the initial multiplicities and matrix digest. This is a certificate for the finite matrix; the required mathematical bridge identifies the literal array count with `v T^(n-1) 1`. The program also checks the finite indices before this certificate implies the published recurrence, the generating-function coefficients through degree 200, the original arrays directly whenever there are at most 16 cells, and Barker's A185761 closed form in `Q(sqrt(17))` for all 200 terms. There are 142 additional applicable cells in the 683-entry parent-table b-file A185769.

These checks support a proposed universal proof. They are not a Lean proof or a claim that finite sequence agreement establishes an unbounded recurrence.

Sources: [A185761](https://oeis.org/A185761), [A185762](https://oeis.org/A185762), [A185763](https://oeis.org/A185763), [A185764](https://oeis.org/A185764), and parent table [A185769](https://oeis.org/A185769), all by R. H. Hardin. `oeis-entries.json` preserves entry revisions and source timestamps. The `b*.txt` files preserve the complete retrieved b-files. OEIS data are attributed under the [OEIS end-user license](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement) (CC BY-SA 4.0); the check implementation is original to this experiment.
