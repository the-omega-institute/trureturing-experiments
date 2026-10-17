# OEIS A400454 and A400825: the test 2^((N-1)/4) ≡ −1 (mod N)

For N = (2^k − 1)·2^(k+j) + 1, the two entries conjecture that 2^((N − 1)/4) ≡ −1 (mod N) implies N is prime. The case j = 3 is A400454 and j = 4 is A400825.

`check.py` (Python 3 with sympy) takes each k with 1 ≤ k ≤ K and checks:
- the implication, with sympy's `isprime` deciding whether N is prime;
- the size bound N < (2^(k+j−1) + 1)², which the proof uses;
- the list of k for which N is prime, against the published terms.

It also reports for how many prime N the hypothesis holds. This bears on the converse, which neither entry conjectures.

Run it with `python3 check.py 400`. Exit 0 and `RESULT PASS` mean the checks hold. `run.log` is the output for K = 400.
