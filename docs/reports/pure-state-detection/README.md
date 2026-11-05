# Pure-state detection under incomplete qubit–qutrit tomography

Source: Veronica Kanachova, Feng Liu and Zhoudunming Tu, *Qubit–Qutrit Quantum
Tomography of hadronic Λϕ and ΛK^{∗0} systems*, arXiv:2609.11151v1,
10 September 2026, Appendix §.7, remark after Proposition 1.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. `output.txt` contains the
result. All computations use integers, rational numbers or exact symbolic
algebra. There is no floating point evidence.

The fidelity check verifies the characteristic polynomial in the restricted
Schmidt basis used in Proposition 1. It reproduces minimum partial-transpose
eigenvalues −6/25 for Schmidt amplitudes (3/5, 4/5) and −1/4 for equal amplitudes.
The remaining quadratic factor has positive trace and determinant.

For arbitrary complex coefficient rows u,v, the symbolic check verifies
η†Mη = −Δ, ψ†η = 0, ‖η‖² = ‖ψ‖², and the three-minor expression for Δ,
where η = (U conjugate(v), −U conjugate(u)) and
Δ = ‖u‖²‖v‖² − |u†v|². It also checks that M is Hermitian.

An independent Gaussian-integer implementation constructs the full 6×6 matrix
using the source's second-factor partial transpose. It checks every six-tuple
in {0,1,−1,i,−i}^6: one zero vector, 744 product vectors and 14,880 entangled
vectors. For each nonzero tuple of squared norm N it verifies the normalized
witness expectation −Δ/N² exactly, including strict negativity in the
entangled cases. Product structure is decided by pivot cross multiplication.

These are checks supporting a preregistration outline. They are not a Lean
proof or a settlement of the published conjecture.
