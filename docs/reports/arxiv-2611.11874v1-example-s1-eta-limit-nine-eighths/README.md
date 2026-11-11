# Corenflos, Example S1: stopping inflation

Source: Adrien Corenflos, [An iterated-logarithm price of anytime validity for log-optimal e-processes](https://arxiv.org/html/2611.11874v1), arXiv:2611.11874v1, Example S1 and equation (4).

Run `python3 check.py` from this directory. The script uses only the Python standard library. It computes the fixed-horizon likelihood ratios for the symmetric endpoint mixture and verifies feasibility against all five nulls for horizons 0 through 200. Given feasibility, the ratio is a numeraire: under the alternative, the expected ratio of any e-variable to it is the mean of that e-variable under the two endpoint nulls.

The script reproduces the exact stopped expectation 30/29 and certifies the reported rounding of the gaps at horizons 50 and 200. It computes each of the five bounded stopping values by exact Snell recursion, checks the proposed inequalities and monotonicity, and cross-checks the first five horizons by enumerating stopping policies on the ordered prefix tree. `check.out` is the recorded output. These finite checks do not prove the all-horizon prediction.
