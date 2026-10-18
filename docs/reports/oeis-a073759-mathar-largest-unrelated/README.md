# OEIS A073759 (Mathar): a(n) = n − lpf(n) for n ≥ 5

Entry: "Largest number < n that is neither a divisor of nor relatively prime to n, or 0 if no such number exists." Conjecture: "a(n) = 2*A046667(n) if n>=5", where A046667(n) = A046666(n)/2 and A046666(n) = n − (smallest prime dividing n).

Route: let p be the least prime factor of n.
- **Composite n ≥ 5.** m = n − p satisfies gcd(m, n) = p < m, so it is unrelated to n.
- **No larger candidate.** A prime dividing both n and some k with n − p < k < n would divide n − k < p, which contradicts the minimality of p.
- **Prime n.** The set is empty, and n − p = 0.

## Program

`check.py` (standard library only). Run `python3 check.py`. It prints `RESULT PASS` when every check holds:
- brute force agrees with 2·A046667(n) for 5 ≤ n ≤ 20000;
- the first 50 OEIS terms agree;
- a(4) = 0 while 4 − lpf(4) = 2.

Reading: `run.log` (exit 0, `RESULT PASS`, 1.4 s).
