# Quantum Volunteer’s Dilemma: fidelity and strategic stability

Source: N. D. Hebdon and D. E. Koh, *Entanglement in the Quantum Volunteer’s Dilemma*, arXiv:2606.08227v2, Eqs. (17), (51), DOI 10.1103/2lb4-m4bg.

Run with Python 3 and SymPy:

```sh
python3 check.py
```

The program uses exact rational and symbolic complex arithmetic. For n = 2, 4, 6 it reproduces Eq. (17) and computes every unilateral payoff from the complete tensor-product circuit with the canonical decoder, original two-parameter strategies, and complementary-bit payoffs. It checks the proposed two-sector identity for every player, including arbitrary diagonal bias and complex coherence. It also checks the scalar optimization algebra, a rational parametrization of the mixed phase-noise family, its partial-transpose characteristic polynomial, and an explicit two-player counterexample above Eq. (51)’s threshold.

The proposed general classification is: on the span of the all-zero and all-one basis vectors, write the resource block as [[a,z],[conj(z),1−a]], put s = (−1)^(n/2) and F = 1/2+s Im(z). For even n ≥ 2, the profile A^n is a Nash equilibrium precisely when z = i s r for some r ≥ 0 and F ≥ 2n/(3n−1). The experiment checks the circuit identity at the listed sizes and the symbolic scalar steps; it is not a Lean proof for arbitrary n.

The scope fixes J† throughout and varies only the prepared density matrix. It makes no assertion about input-adapted decoders, retuned strategy profiles, or approximate Nash equilibria. The source’s matched-decoder pure-family theorems are preserved.

`output.txt` contains the successful output, ending in `ALL_CHECKS_OK`.
