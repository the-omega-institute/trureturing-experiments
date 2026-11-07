# Virasoro singular vector and VBKM Cartan cancellation

Source: Vladyslav Zveryk, *Virasoro extensions of generalized Kac-Moody algebras*, arXiv:2609.15729v2, Conjecture 4.8, Definitions 2.1 and 3.1, and equation (3.1).

Run `python3 check.py`. Only Python's standard library is required. The output is stored in `check.out`; a successful run ends with `ALL_CHECKS_OK`.

The fidelity checks evaluate the source's singular-weight formula at central charge 4 and indices (2,2), obtaining -3/8, and reproduce its normalized bracket [L₂,L₋₂]=4B+6C, with conventional charge c=12C. The independent rational Verma-module calculation checks the level-four word

P = L₋₁⁴ − 3L₋₂L₋₁² − (7/4)L₋₂² + 3L₋₃L₋₁ + (3/2)L₋₄.

Positive modes 1 through 8 annihilate P applied to the highest-weight vector. The joint L₁/L₂ matrix at level four has rank four on five PBW monomials. There are 300 consistency checks of the Virasoro commutator on basis words of levels 0 through 4 and modes -2 through 2. The implementation computes required intermediate words rather than setting words beyond a cutoff to zero.

For q=1 through 100, the program verifies the proposed Cartan matrix, determinant, central and degree pairings, condition (4.3), a strictly positive root-separating functional, and negative rank-one highest weights. These are finite checks of the family formulas. The iterated adjoint action of P on H₁ is zero: its last two contributions are −6δL₋₄ and +6δL₋₄, with δ=−1/3.

A separate negative-Witt action L₋ᵣ(tᵏ)=−ktᵏ⁻ʳ gives P(t)=(1/4)t⁻³. This detects a nonzero negative word without assuming PBW independence.

The recursive representation of the full q=1 presentation is also checked on negative words of length at most two, modes -3 through 3, and two Cartan highest weights. All 3,654 exact generator-relation checks pass. The expanded adjoint singular vector is annihilated by both raising generators and positive modes 1 through 5 in these representations. These finite checks support the recursive construction; they do not replace its universal verification.

This program does not prove triangular decomposition, nonzero survival in the defining quotient, the universal family statement, or ideal/Cartan disjointness. Those require the presentation and weight arguments in the preregistration. No Lean settlement proof is included.
