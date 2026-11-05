# Exact checks for the proposed F_sigma classification

`check.py` reproduces the source norm identity for the finite-rank operator
`T_{3,2}` at `p = 2`: the exact value satisfies
`||T_{3,2}||^2 = 3/4`, using the binomial torus moments
`E|z_1+z_2|^(2k) = binomial(2k,k)`. It also checks the proposed dense
nonclosed example `B = Q ∩ [1,∞)` against the necessary upward-closure
property of one common operator: `2 ∈ B` while `2 < sqrt(5)` and
`sqrt(5) ∉ B`. All computations use integers, `Fraction`, and exact SymPy
objects; no floating-point value is used as evidence.

The captured output is in `output.txt`.
