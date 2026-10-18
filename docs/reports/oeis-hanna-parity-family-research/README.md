# Hanna OEIS parity and integrality conjectures: A393173, A393850, A382122, A378258, A375449, A375452, A365574, A380711, A379764, A380681–A380686

These entries' conjectures are settled by a GPT Pro research answer. The orchestrator checked every step, and this program recomputes each sequence from its literal definition (or the exact rewrite of a bilateral sum) and tests the conjecture.

## Program

`check.py` (standard library only; reads the saved entry texts). Run `python3 check.py`; it prints `RESULT PASS` when every check below holds.
- **Each entry.** The first 12 terms match the OEIS data, all computed coefficients are integers, and the conjecture holds on the stated range:
  - n ≤ 40 for A393173, A393850 (n ≥ 1), A382122, A378258, A380711 and A379764 (both conjectures);
  - 2 ≤ n ≤ 60 for A375449 and A375452;
  - 0 ≤ k ≤ 59 for A365574 (Conjecture 1).
- **A393850.** a(0) = 1 is odd and 0 is not a power of 2, so the literal comment fails only at n = 0.
- **A382122.** A second solution branch of the absolute-value equation is integral with the same parity.
- **A365574.** The explicit reversion coefficients r_k = ((2 − k) + (k + 1)(−2)^{k−1})/3 satisfy the definition for n = 2, …, 19.
- **A380681–A380686.** The bilateral equation holds to order 9–12.
  - For p = 2, 3, 5, 7, 11, 13 the solution is integral and matches the entries.
  - For p = 4, 6, 8, 9, 10, 15 the first non-integer coefficient is at index q − 1 (q the least prime factor) and equals Π_{j<q}(1 + jp)/q!.
- **The bilateral identity.** Ψ₁(x, Y) = 0, and Ψ₃(x, Y) ≡ 0 (mod 3), for random integer Y.

Reading: `run.log` (exit 0, `RESULT PASS`, 10.8 s).
