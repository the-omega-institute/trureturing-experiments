# OEIS A361032: integrality of Bala's quartic super-Catalan table

A361032 is the table T(n,k) = F(n)·(4k)!/(k!·(k+n+1)!³) with F(n) = (4n+4)!/(8·(n+1)!). The entry says that this choice of F(n) "appears to produce integer values".

`check.py` (Python standard library) checks, for n, k ≤ M:
- every T(n,k) is an integer;
- the 2-adic valuation used in the proof holds: v₂(T(n,k)) = 3·(s₂(n+k+1) − 1);
- the floor inequality ⌊4a⌋ + ⌊4b⌋ ≥ 3⌊a + b⌋ holds on rationals, which is the odd-prime step;
- the published antidiagonals match.

Run it with `python3 check.py 45`. Exit 0 and `RESULT PASS` mean every check holds. `run.log` is the output.
