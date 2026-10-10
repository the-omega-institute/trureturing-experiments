# Krishna's C*-algebraic Novak conjecture: exact counterexample in M₂(ℂ)

`check.py` (SymPy, exact) works in A = M₂(ℂ) with n = 3 and the same self-adjoint coordinate in every position l:
- A₁ = 3πI, A₂ = π·diag(5, 1), A₃ = (π/4)[[19, √15], [√15, 5]];
- cos x := (exp(ix) + exp(−ix))/2.

It checks that (A_j − A_k)² is 4π²I, 4π²I and π²I, so the cosines are I, I and −I. It then checks that, for d = 2, 3 and 5, the matrix N = [∏_l (1 + cos(x_{j,l} − x_{k,l}))/2 − 1/n] has scalar block (1/3)[[2,2,2],[2,2,−1],[2,−1,2]] and ⟨Nx, x⟩ = −2·I for x = (−2, 1, 1)·1_A. `check.log` holds the run (exit 0).
