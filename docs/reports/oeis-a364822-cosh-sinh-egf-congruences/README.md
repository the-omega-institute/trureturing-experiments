# OEIS A364822: Tebni's congruences for cosh(x)/(1 − 2 sinh(x))

The conjectures checked are:
- **C1.** a(p) ≡ 2 (mod p) for odd primes p.
- **C2.** a(2^m) ≡ 1 (mod 2^m).

C2 is checked in the stronger form v₂(a(n) − 1) = v₂(n) + 2 for every even n ≥ 2.

The proofs use two expansions of the e.g.f.:
- **C1.** a(n) = Σ_{m≤n} Σ_j C(m,j)(−1)^j((m−2j+1)^n + (m−2j−1)^n)/2. By Fermat this reduces modulo p to Σ_m Σ_j C(m,j)(−1)^j(m−2j) = 2.
- **C2.** a(n) = Σ_{r odd} (1/(2r)) Σ_j C(r,j)(−1)^j(r−2j)^(n+1).
  - Here the r = 3 term is 3^n − 1, with v₂ = v₂(n) + 2 by the lifting-the-exponent lemma.
  - For odd r ≥ 5 the terms are ≡ 0 modulo 2^(v₂(n)+3). This rests on a parity lemma: the number of j < r/2 with C(r,j) odd and r − 2j ≡ ±3 (mod 8) is even.

`check.py` (Python 3 with sympy) checks:
- both expansions against a series expansion and the published data;
- C1 for odd p ≤ NMAX;
- the 2-adic statement for even n ≤ NMAX;
- the parity lemma for r ≤ R.

Run `python3 check.py 256 100001`. Exit 0 and `RESULT PASS` mean every check holds (about 35 s).
