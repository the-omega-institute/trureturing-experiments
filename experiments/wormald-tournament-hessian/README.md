# Wormald, Tournaments with many Hamilton cycles (arXiv:2610.16245), Conjecture 6

Exact rational computation of 𝓗 (determinant of the Hessian at θ = 0 of θ ↦ log h(z e^{iθ})) from the paper's definition, by implicit differentiation of K(y) = det(I − (½I + A) diag(y)) at y = (2/d, …, 2/d); K's partial derivatives are exact determinants (multilinearity in the columns).

    python3 wormald-hessian.py       # standard library: reproduces the paper's values √(d/𝓗) = 2, 44/19, 3032/1231 for the rotational tournaments R_3, R_5, R_7
    python3 wormald-closed-form.py   # needs sympy: random regular tournaments d = 3..11; checks d/𝓗 = (det(dI+S) / (det(S+idI)/i))^2 with S = A − Aᵀ

`wormald_hessian_core.py` is the shared module (identical to `wormald-hessian.py`).
