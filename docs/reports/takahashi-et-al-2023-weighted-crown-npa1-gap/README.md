# Takahashi et al.: the first-level projector SDP is not exact on the open weighted-crown interval

`check.py` supports a proof for the question left open in Appendix C.1 of arXiv:2307.15688v3 (Takahashi, Rayudu, Zhou, King, Thompson, Parekh). The question: is NPA₁(Proj) inexact for the weighted crown H = x h_ab + Σᵢ (h_ai + h_bi) (vertices a, b, 1..n, n ≥ 3) on the interval L < x < U, with L = (n+2)²/(4(n+1)) and U = (n+2)/3?

The checks support the answer yes. For every such n and x, an explicit real moment matrix is feasible, with objective n + 1 + (x − L)²/(32n) > n + 1 = λ_max(H).

- (1) For n = 3..10, at five rational points of (L, U) each, the moment matrix over {I} and all vertex-pair variables satisfies the level-one Proj constraints: M(I,I) = 1, idempotence, symmetry, and the anticommutation relation for every path i–j–k. It is positive semidefinite by exact rational elimination. For n = 3: 7 positive pivots and 4 zero directions. Its objective equals n + 1 + (x − L)²/(32n).
- (2) Symbolic in (n, x): the objective identity, the Schur complement entry 2n p(p₀ − p), and det S = 2n q³(4p₀ − c(c − 1/(2n))). A numerical spectrum check covers n = 5.
- (3) Exact diagonalization for n = 3, 4, 5 at three points of (L, U) gives λ_max(H) = n + 1, as in the source's Eq. (C2).
- (4) The interval is nonempty exactly for n ≥ 3.

Run: `python3 check.py` with numpy 2.4 and sympy 1.14 (exit 0; final line `bad= 0`; about 5 s).
