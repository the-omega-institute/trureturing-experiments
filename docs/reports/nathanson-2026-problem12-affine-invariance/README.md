# Nathanson (arXiv:2604.04781), Problem 12: H(A_q) = H(sA_q + t) for every nonzero s

Source: "Let (A_q) and (A′_q) be asymptotically strictly decreasing sequences of sets of integers such that A′_q = sA_q + t for all q ∈ ℕ and some s, t ∈ ℤ with s ∉ {0, ±1}. Does this imply H(A_q) = H(A′_q)?"

Route: the answer is yes, for every nonzero integer s.
- **Injectivity.** f_h(x) = sx + ht is injective on ℤ, and injective maps commute with intersections over a nonempty index set.
- **Sumsets.** h(sB + t) = f_h(hB).
- **Conclusion.** hA′ = f_h(hA) and ⋂ hA′_q = f_h(⋂ hA_q), so the two equalities are equivalent.

The paper's Theorem 15 covers only units ε = ±1. The proof above needs injectivity, not invertibility.

## Program

`check.py` (standard library only). Run `python3 check.py`; it prints `RESULT PASS` when every check holds.
- **The check.** It compares H(A_q) with H(sA_q + t), for 1 ≤ h ≤ 4, on random finite families of finite integer sets. The multipliers are s ∈ {−3, −2, 2, 3, 5} and the shifts t ∈ {−4, 0, 7}.
- **The control.** It reports how many families have nontrivial H, so the check is not vacuous.

Reading: `run.log` (exit 0, `RESULT PASS`). There are 6000 cases and 0 mismatches; 291 families have H missing some h ≤ 4.
