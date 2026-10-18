# Maity–Furuichi–Dutta (arXiv:2607.01805): the ordering B₁ ≤ B₂ fails for non-commuting density operators

Source: Eq. (9), "We observe that for non-commutative density operators ρ and σ, D_α(ρ||σ) ≤ B₁(ρ||σ) ≤ B₂(ρ||σ)". Section 4, Case 5 repeats it: "for any two operators ρ and σ, B₁(ρ||σ) is a tighter bound of Tsallis relative entropy than B₂(ρ||σ)".

**Counterexample.** α = r = 1/2, z = 2, q = 1, ρ = diag(4/5, 1/5) and σ = (1/125)[[28, 6√6], [6√6, 97]]. Both are density matrices with spectrum {4/5, 1/5}, and they do not commute.
- B₁ = 4(1 − √((2117 + 288√2)/3125)) ≈ 0.404950.
- B₂ = 2 − 20/√154 ≈ 0.388354.
- D_α = 48/125 = 0.384.

So 48/125 < B₂ < 2/5 < B₁: both are upper bounds, but B₂ is the tighter one.

**Mechanism.** The paper's Theorem 2 gives B₂ ≤ B₁ for commuting pairs, strictly when ρ ≠ σ. A small rotation of σ keeps the strict inequality, so non-commuting pairs near any commuting pair with ρ ≠ σ violate B₁ ≤ B₂.

## Program

`check.py` needs sympy and mpmath (`python3 -m pip install sympy mpmath`). Run `python3 check.py`; it prints `RESULT PASS` when every check holds:
- it evaluates the counterexample at 50 digits and checks it against the closed forms;
- it checks B₁ > 2/5 > B₂ > D_α = 48/125;
- it checks B₁ > B₂ on six near-commuting pairs (rotations by 0.01 and 0.1).

It also reports how closely the same code reproduces the paper's four Case 3 numbers; this is reported only, not a pass condition. Three agree to within 5·10⁻⁶. The fourth, B₁ at α = 0.99, is 0.252973 here against the paper's 0.252861.

Reading: `run.log` (exit 0, `RESULT PASS`).
