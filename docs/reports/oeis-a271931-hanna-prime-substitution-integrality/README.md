# OEIS A271931 / A271932 / A271933 (Hanna 2016): integer solutions of G(x) = G(x^p + p x G(x)^p)^(1/p) for prime p

Conjecture (in all three entries): "for prime p there exists an integer series G(x) that satisfies: G(x) = G( x^p + p*x*G(x)^p )^(1/p) with G(0)=0, G'(0)=1".

## Route

- **Reversion.** R = G⁻¹ satisfies R(y^p) = R(y)^p + p·y^p·R(y).
- **Recursion.** The coefficient of y^{n+p−1} gives p·r_n = [y^{n+p−1}](R(y^p) − R_{<n}(y)^p − p·y^p·R_{<n}(y)).
- **Divisibility.** The bracket is divisible by p, because R_{<n}(y)^p ≡ R_{<n}(y^p) (mod p) and (n + p − 1)/p < n.
- **Conclusion.** R and G are integral.

## Program

`check.py` (standard library only; reads the three OEIS entry texts saved here). Run `python3 check.py`; it prints `RESULT PASS` when every check below holds.
- **Primes p ≤ 23.** Every division by p is exact, and the inverse G satisfies G^p = G(x^p + p·x·G^p) to the computed order (30 for p ≤ 5, down to 12 for p = 23).
- **The OEIS entries.** For p = 5, 7, 11, G agrees with the first 12 terms of A271931, A271932 and A271933.
- **Composite exponents.** For p = 4, 6, 8, 9 the reversion has a non-integer coefficient (r₃ = −1/2, −3/2, −5/2; r₄ = −35/3).
- **Catalan.** For p = 2 and 3, G is the Catalan generating function, as the entries' comparison states.

Reading: `run.log` (exit 0, `RESULT PASS`).
