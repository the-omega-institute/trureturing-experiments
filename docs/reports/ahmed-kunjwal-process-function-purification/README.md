# Ahmed–Kunjwal process-function purification check

This report checks the finite reversible extension and purification construction from Ahmed and Kunjwal, *Characterizing unitaries via quasi-process functions*, arXiv:2610.00579v1, Appendix A, against Conjecture VII.1. `check.py` uses only exact integer, rational, and SymPy exact complex arithmetic.

The script exhaustively enumerates one and two party binary/singleton process functions, checks the published three-party AF/BW process and four causal three-party examples, verifies the reversible map and its inverse, verifies the source preparation/future-discard reduction, tests the multilinear normalization identity on exact complex matrices with constant column sums, and checks linked Kraus completeness for identity and dephasing local CPTP families. The finite checks are evidence for the stated route and do not settle the universal conjecture.

Run:

```sh
python3 check.py
```

The saved output is `check.out`; its final line is `ALL_CHECKS_OK`.
