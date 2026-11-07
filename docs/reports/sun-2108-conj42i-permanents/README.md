# Sun's absolute-difference permanent congruences

Source: Zhi-Wei Sun, *Arithmetic properties of some permanents*, arXiv:2108.07723v7, Conjecture 4.2(i), equations (4.3) and (4.4), and Remark 4.2. The source was first posted in 2021; v7 is dated 6 June 2022.

Run with Python 3, without external dependencies:

```sh
python3 docs/reports/sun-2108-conj42i-permanents/check.py
```

`output.txt` contains the complete deterministic output. The program exits 0 and ends with `ALL_CHECKS_OK`.

- Fidelity: the two q-integer determinant identities stated in Remark 4.2 hold for orders 2 through 9 and q = 0, 1, 2, 3. Bareiss elimination uses exact integer division and is checked against the Leibniz formula at orders 2 through 6.
- Target: direct exact subset dynamic programming checks both permanent congruences for every odd prime through 19. Direct permutation sums independently agree through 7; Ryser's inclusion-exclusion independently agrees through 13.
- Reduction: 151986 entries verify the alternating row and column sign identity for both matrices at every odd prime through 101. The alternating row sums also have the asserted exact integer values.
- Classical ingredient: the circulant permanent congruence is checked for every function from F_p to F_p at p = 3, 5 (3152 functions in total). Direct enumeration confirms that the permutations commuting with the one-step cycle are exactly the translations at p = 3, 5, 7.

All evidence uses integers, including integer modular arithmetic. These are finite checks supporting a preregistration draft. They establish no universally quantified theorem and contain no Lean proof. Conjecture 4.2(ii) is outside the scope.

Producer: codex-cli research flight `b588670809297ce87877c2aa`, with prior disclosure `repo-prior-exposed`. The independent implementation uses the candidate's disclosed proposed mathematical route.
