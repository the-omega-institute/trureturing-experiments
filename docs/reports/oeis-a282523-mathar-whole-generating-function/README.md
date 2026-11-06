# Mathar's OEIS A282523 generating function

The source is [OEIS A282523, revision 9](https://oeis.org/A282523/internal). R. H. Hardin defines the sequence of labelled n by 3 binary arrays whose occupied cells have at most two occupied king-move neighbors; the offset is 1. R. J. Mathar's FORMULA line conjectures the complete rational generating function. The parent is [A282528](https://oeis.org/A282528/internal). OEIS source data are available under [CC BY-SA 4.0](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).

Run `python3 check.py` with Python 3.10 or later. It uses only the standard library and exact integers. `output.txt` is its stdout; successful execution ends with `ALL_CHECKS_OK`.

The program independently implements the literal cell-neighborhood definition and the middle-row rule, checks every one of 512 row triples, verifies the proposed 64-to-13 state quotient, and checks the complete polynomial certificate coefficientwise. It also verifies the degree-eight annihilator directly on the raw 64-state transfer, reproduces all 21 displayed source terms, exhaustively enumerates every binary board for heights 1 through 6, and compares raw and quotient transfers through height 210. The proposed certificate is attributed to search seat codex-search-b; this program supplies an independent implementation of its checks.

The intended proof must establish the bijection from literal arrays to accepted row words for every positive height and propagate the finite certificate by induction/formal-series uniqueness. These finite checks do not establish that bridge in Lean. No Lean proof of the conjecture or its negation is included. Generic transfer-matrix rationality alone does not determine the conjectured numerator and denominator.
