# OEIS A207969/A207970 and A208034 (Hanna 2012): exp(Σ 5F_n^{2k}xⁿ/n) and exp(Σ 2P_n^{2k}xⁿ/n) are integer series

Conjectures:
- **A207969, A207970.** "exp( Sum_{n>=1} 5*Fibonacci(n)^(2*k) * x^n/n ) is an integer series for integers k>=0."
- **A208034.** For all k ≥ 1: (1) exp(Σ 2·Pell(n)^{2k}xⁿ/n) is an integer series; (2) with 2·Pell(n)^{2k−1} it is not; (3) with Pell(n)^{2k} it is not.

## Route

exp(Σ c_n xⁿ/n) is integral iff c_{pm} ≡ c_m (mod p^{v_p(m)+1}) for all primes p and all m ≥ 1 (Dwork). Write the even powers as integer combinations of traces of powers of integer matrices:
- 5^k·F_n^{2k} = (φⁿ − φ′ⁿ)^{2k};
- 8^k·P_n^{2k} = (αⁿ − βⁿ)^{2k}.

Such combinations satisfy the congruences. Dividing by 5^{k−1} (respectively 2^{3k−1}) matters only at 5 (respectively 2). There both sides of every congruence are divisible by the needed power, because v₅(F_n) = v₅(n) and v₂(P_n) = v₂(n).

The failures are explicit:
- 5F^{2k+1} fails at p = 3 (c₃ − c₁ ≡ 2 (mod 3));
- F^{2k} with k ≥ 1 fails at p = 5;
- 2P^{2k−1} fails at p = 3;
- P^{2k} fails at p = 2.

## Program

`check.py` (standard library only; reads the saved entry texts). Run `python3 check.py`; it prints `RESULT PASS` when every check below holds.
- **Even powers.** exp(Σ 5F^{2k}xⁿ/n) for k ≤ 8 and exp(Σ 2P^{2k}xⁿ/n) for 1 ≤ k ≤ 8 are integral to order 80. Their local congruences hold for p ≤ 13 and m ≤ 40.
- **The failing cases.** Each one is non-integral, with the stated witness.
- **Valuations.** v₅(F_n) = v₅(n) and v₂(P_n) = v₂(n) for n ≤ 160.
- **The entries.** The first 12 terms of A207969, A207970 and A208034 match.

Reading: `run.log` (exit 0, `RESULT PASS`).
