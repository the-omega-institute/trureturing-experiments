# OEIS A394815: twin-prime centers with gcd(k, σ(k)) = 2 are divisible by 36

A394815 lists the averages k of twin prime pairs with gcd(k, σ(k)) = 2. Its entry conjectures that every term is divisible by 36.

`check.py` (Python 3 with sympy) runs over every twin-prime center k < N and checks:
- the conjecture;
- the published terms;
- the proof's 2-adic step: if v₂(k) = 1 and 3 | k, then 3 | σ(k).

Run it with `python3 check.py 2000000`. Exit 0 and `RESULT PASS` mean every check holds. `run.log` is the output for N = 2·10⁶ (91 terms).
