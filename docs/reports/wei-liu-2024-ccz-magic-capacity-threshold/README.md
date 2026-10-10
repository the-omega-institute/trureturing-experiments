# Wei–Liu: the CCZ magic capacity threshold under local depolarizing noise is not 1/3

`check.py` checks a refutation of the conjecture in Section VII of arXiv:2410.21215v1: "We conjecture that the magic capacity threshold for CCZ under local depolarizing noise is 1/3."

1. It enumerates all 1080 pure three-qubit stabilizer states by breadth-first search over the Clifford generators H, S and CNOT, and checks that the maximum squared overlap with ψ = CCZ|+++⟩ is 9/16.
2. For the six-qubit stabilizer input Ω = 8^{-1/2} Σ_x |x⟩_A|x⟩_B and the witness W = V†[(9/16 I − |ψ⟩⟨ψ|) ⊗ |000⟩⟨000|]V, with V the CNOT layer A_i → B_i, it checks two things:
   - Tr(W ρ_λ) = (9/16)(1 − λ/2)³ − (1 − 3λ/4)³ against an exact formula. At λ = 1/2 the value is −7/1024.
   - The formula is negative at every λ = k/1000 in [0, 1/2], in exact rationals.
3. It checks ⟨φ|W|φ⟩ ≥ 0 on 3000 random six-qubit stabilizer states.
4. It checks that the noise can act before or after CCZ without changing the witness values.

Run: `python3 check.py` (exit 0, final line `ALL_OK_BOTH_ORDERS`).
