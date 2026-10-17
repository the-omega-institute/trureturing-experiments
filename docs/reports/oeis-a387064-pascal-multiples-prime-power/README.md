# OEIS A387064: multiples of p^k in rows 0..p^k of Pascal's triangle

A387064(n) is the number of entries in rows 0..n of Pascal's triangle that are multiples of n. The entry conjectures a(p^k) = (p − 1)·p^(k−1) for primes p.

`check.py` (Python 3 with sympy) checks:
- the published terms for n ≤ 25;
- the conjecture for every prime power p^k ≤ N;
- the two facts the proof uses: no entry of a row m < p^k is divisible by p^k, and C(p^k, j) is divisible by p^k exactly when p ∤ j.

Run it with `python3 check.py 300`. Exit 0 and `RESULT PASS` mean every check holds. `run.log` is the output for N = 300.
