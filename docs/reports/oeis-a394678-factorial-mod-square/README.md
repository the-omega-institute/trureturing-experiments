# OEIS A394678: when (n! − 1) mod n² is prime

A394678 is a(n) = (n! − 1) mod n². Its entry conjectures that, apart from a(4) = 7, a(n) is prime only when n is a prime with n² − n − 1 prime.

`check.py` (Python 3 with sympy) takes each n ≤ N and checks:
- the conjecture;
- the published terms;
- the two facts the proof uses: a(p) = p² − p − 1 for prime p (Wilson), and n² | n! for composite n > 4, so that a(n) = n² − 1.

Run it with `python3 check.py 3000`. Exit 0 and `RESULT PASS` mean every check holds. `run.log` is the output for N = 3000.
