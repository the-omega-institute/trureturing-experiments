# OEIS A361033, A361034, A361035: parity of Bala's quartic super-Catalan analogues

Peter Bala's sequences are a_j(n) = C_j·(4n)!/(n!·(n+j)!³) with C_1 = 3, C_2 = 2520 and C_3 = 9979200. He conjectures that a_j(n) is odd iff n + j is a power of 2.

`check.py` (Python standard library) checks, for n ≤ N and j = 1, 2, 3:
- the quotients are integers;
- the parity conjectures hold;
- the published data match;
- the 2-adic valuation used in the proof holds: v₂(a_j(n)) = 3·(s₂(n + j) − 1), where s₂ is the binary digit sum.

Run it with `python3 check.py 600`. Exit 0 and `RESULT PASS` mean every check holds. `run.log` is the output.
