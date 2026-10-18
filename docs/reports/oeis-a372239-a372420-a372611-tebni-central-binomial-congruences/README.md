# OEIS A372239, A372420, A372611: Tebni's congruences at primes

Mélika Tebni's sequences are a_c(n) = [x^n] (1 + c·x)/((1 − 2x)·√(1 − 4x)) for c = 2, 1, 3. Each entry makes two conjectures:
- for Pythagorean primes p, a_c(p) ≡ 4 + c (mod p);
- for primes p ≡ 3 (mod 4), a_c(p) ≡ −c (mod p).

`check.py` (Python 3 with sympy) checks:
- the published data;
- all six conjectures for primes 3 ≤ p < P;
- the proof's two congruences b(p) ≡ 2((−1)^((p−1)/2) + 1) and b(p−1) ≡ (−1)^((p−1)/2) (mod p), where b(n) = Σ_{k≤n} C(2k,k)·2^(n−k).

Run it with `python3 check.py 1200`. Exit 0 and `RESULT PASS` mean every check holds. `run.log` is the output.
