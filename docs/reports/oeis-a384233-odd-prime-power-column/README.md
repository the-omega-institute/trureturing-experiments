# OEIS A384233: column 2 is the odd prime powers

A384233 is the square array whose entry T(n,k) is the n-th number whose largest odd noncomposite divisor (1 or an odd prime) is its k-th divisor. Its entry conjectures T(n,2) = A061345(n): column 2 consists of the powers of odd primes, in increasing order.

`check.py` (Python 3 with sympy) checks, for every m < M:
- position 2 occurs exactly for the odd prime powers;
- the array, rebuilt from these positions, matches the published antidiagonals;
- column 2 matches A061345, which has offset 0 (1, 3, 5, …).

Run it with `python3 check.py 20000`. Exit 0 and `RESULT PASS` mean every check holds. `run.log` is the output.
