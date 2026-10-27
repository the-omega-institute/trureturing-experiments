# Hardin A204612 recurrence checks

Run `python3 check.py` in this directory with Python 3. The program uses integers only and needs no third-party package or network connection. Its complete output is in `output.txt`; success ends with `ALL_CHECKS_OK`.

The source is R. H. Hardin's OEIS A204612, submitted 17 January 2012, revision #9 of 3 October 2025. The literal array definition and empirical recurrence are at <https://oeis.org/A204612/internal>. The linked table is <https://oeis.org/A204612/b204612.txt>. `source-data.json` records the definition, formula, 16 displayed values and all 210 linked values used as fidelity observations. Source-data reuse is subject to the OEIS End-User License Agreement, <https://oeis.org/LICENSE>.

The program constructs all terminal states and transitions from actual pairs of binary rows. A separate representation enumerates every two-, three- and four-row binary matrix directly. It compares the transfer counts with all source data, verifies the specified nonzero polynomial residual, checks the images of every generator of its invisible invariant subspace, and checks the recurrence for every index 16 through 500. It also checks an auxiliary degree-four vector identity.

These are finite exact checks of the proposed all-height argument. They do not substitute for the append/delete bijection and the induction propagating the invariant to every height, and no Lean proof is included. The recurrence is never used to define the counted objects or generate their counts.
