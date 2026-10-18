# OEIS A344695: gcd(sigma(n), psi(n))

Check for the open-problem preregistration on the Karttunen–Munn conjecture in OEIS A344695: the primes that occur as values of gcd(σ(n), ψ(n)) are exactly the odd primes.

Run: `python3 check.py 2000000 200` (Python 3, standard library only).

It checks:
- the published data;
- that a(n) ≠ 2 for n ≤ 2·10^6, together with the dichotomy used in the proof: a(n) is odd for square n, and (p + 1) | a(n) for a prime p of odd exponent otherwise;
- the construction for every odd prime q ≤ 200. Take r the least prime ≡ 1 (mod q(q+1)(1+q+q²)) with r > 1 + q + q², and n = q²·r^(q−1); then gcd(σ(n), ψ(n)) = q.

It also lists the first n with a(n) = q for the primes q ≤ 23 below the sieve bound.

`run.log` holds the output of the run cited by the preregistration.
