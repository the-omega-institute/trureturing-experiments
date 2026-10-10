# Aloy–Müller-Rigat–Tura–Fadel: classical validity of the five proposed three-outcome PI Bell inequalities

`check.py` checks the five 3PIBIs proposed "for any N > 3" in Table III of arXiv:2406.11792v1 (Entropy 26 (2024) 816, Table 3): B = α₁P̃₀ + α₂P̃₀₀ + α₃P̃₀₁ + α₄P̃₁₀ + α₅P̃₁₁ + β_c ≥ 0 with rows (1,1,0,−2,0,0), (1,1,−2,−2,2,0), (−2,1,2,2,0,4), (−6,1,4,4,2,12), (−6,1,4,0,0,24).

- For 1 ≤ N ≤ NMAX it evaluates each row on every count vector c_{a,a'} (Table II of the source) and checks that the minimum is 0 for every N ≥ 2 (nonnegative for N = 1).
- For N ≤ 4 it enumerates every local deterministic strategy and recomputes the observables from their definitions P_{a|x} = Σ_i p(a_i|x_i) and P_{ab|xy} = Σ_{i≠j} p(a_i|x_i)p(b_j|y_j), checking agreement with the count formulas.
- It checks symbolically (sympy) the identities of a proof route for rows 2–5: B₂ = (R−S)² + 2K, the identity for 2B₃, B₄ = G(K,u,v) and B₅ = G(c₀₁,·,·) + G(c₁₀,·,·) with G(k,u,v) = 6(k−1)(k−2) + w² + (6k−7)w + 2uv, w = u+v.

Run: `python3 check.py 12` with sympy (exit 0; final line `NMAX=12 minima(N=1..4)=[[0, 0, 0, 0, 12], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]] identities_ok=True bad=0`).
