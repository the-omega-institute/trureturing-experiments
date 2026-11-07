# OEIS A240515 order-48 recurrence check

This report independently checks R. H. Hardin, OEIS A240515, revision 10 (23 July 2025), whose name is “Number of nX4 0..1 arrays with no element equal to exactly two horizontal and vertical neighbors, with new values 0..1 introduced in row major order.” The parent table is A240519, whose `k=4` entry is `[order 48]`. The source has no arXiv identifier.

`check.py` uses only exact Python integers. It exhaustively counts the four smallest open-boundary boards, reproduces all 49 values needed by the recurrence boundary from an independent 256-state row-pair transfer, checks the published recurrence for every available b-file value n=49..210, and verifies the order-48 annihilator on every row-pair coordinate. The transfer has 256 states and 626 directed transitions; all 256 residuals are zero. The embedded table is the A240515 b-file downloaded from OEIS during this preregistration.

Run:

```sh
python3 check.py
```
