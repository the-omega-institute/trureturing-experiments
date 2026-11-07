# Positive definite symmetric functions without a common Reznick exponent

Source: Sebastian Debus and Robin Schabert, *Any-dimensional Positivstellensätze for symmetric functions*, arXiv:2606.27551v2, Remark 6.7 and Section 8; the fidelity identity is from Example 6.6.

Run `python3 check.py` with Python 3 and SymPy. The stored standard output is `check.out`; successful completion ends with `ALL_CHECKS_OK`. All mathematical comparisons use integer or rational arithmetic.

The program verifies the source's symbolic deletion identity for its decic A, the Robinson nonnegativity decomposition, and the identity between the proposed degree-18 symmetric function H and the sum of weighted Robinson triples in dimensions 1 through 7. It constructs the complete 55-monomial ternary degree-9 moment matrix from 12 rational moments, checks its parity decomposition, and verifies exact positive-pivot LDLᵀ reconstructions for the two representative blocks of sizes 15 and 10. The other two blocks are congruent copies.

The functional gives L(H) = −3231/10⁹ and L(p₁₈) = 3. Thus it separates every H + εp₁₈ with 0 < ε < 1077/10⁹ from the ternary SOS cone. The fixed witness ε = 10⁻⁷ has value −2931/10⁹. Exact small substitution and coefficient checks also support the proposed removal of a multiplier uniform in dimension.

The general argument uses positivity of Robinson triples, convergence of ordinary power sums on ℓ², closedness of each fixed-variable homogeneous SOS cone, and extraction of the highest coefficient in an appended variable. These infinite and general statements are mathematical proof obligations for the subsequent formalization; this program does not prove them in Lean. In particular, finite checks alone do not establish failure of every common exponent.

The two representative LDL certificates contain 350 rational entries including zeros, at most 261 numerator bits and 274 denominator bits, with 10,971 bytes in their rational string representations. These are arithmetic size measurements; Lean elaboration and kernel checking costs have not been measured.
