# Gnedin predictive-rule splitting

`check.py` performs exact algebra checks for Alexander Gnedin, “Exchangeable
Testing Against an Unknown Benchmark”, arXiv:2608.04838v2 (25 October 2026),
Section 6.1. It reproduces the signed cubic factorisation from the paper and
checks the triangular-grid zero and largest-moved-index calculations used by
the proposed all-degree argument. SymPy expressions are over exact rationals;
the permutation checks cover every permutation for degrees `d = 2, ..., 9`.

`check.out` is the captured stdout from Python 3 running `check.py`; its final
line is `ALL_CHECKS_OK`.
