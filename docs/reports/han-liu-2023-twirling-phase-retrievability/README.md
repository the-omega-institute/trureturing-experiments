# Han–Liu: the twirling phase-retrievability formula fails with multiplicities

`check.py` supports a refutation of the final Conjecture of Section 5 of arXiv:2311.16349v1 (Rep. Math. Phys. 93 (2024) 87–102). The conjecture: for π = m₁π₁ ⊕ … ⊕ m_dπ_d with pairwise inequivalent irreducible πᵢ and I_k ≤ d < I_{k+1}, pr(Φ_π) = max{β(Φ_π), k}. Here β(Φ_π) = max mᵢ (the source's theorem) and I_k is the minimal length of a phase-retrievable frame for ℂ^k.

The counterexample family: the three characters of C₃, each with multiplicity 2s. Then d = 3, k = 1 and β = 2s, but pr ≥ 3s. More precisely, for a three-block pinching with blocks m₁, m₂, m₃ (N = Σmᵢ, r = max mᵢ), pr = max{r, ⌊N/2⌋}.

- (1) The C₃ twirl equals the three-block pinching. Checked numerically on random complex matrices for four block triples.
- (2) The encoding V(u,v,w) = ((u,v),(u,w),(v,w))/√2 is an isometry, and X ↦ Φ(VXV*) is injective on all operators (rank n² over ℂ). Checked for all 24 block triples with N ≤ 12 and ⌊N/2⌋ > max mᵢ, including (2s,2s,2s) for s = 1, 2, 3.
- (3) The upper-bound mechanism: for random subspaces of dimension max{r, ⌊N/2⌋} + 1, a vector x of M ∩ D⁻¹M and Dx collide under the pinching, with D = diag(1, ω, ω²) on the three blocks.
- (4) I₂ ≥ 4: for three vectors in ℂ² the intensity kernel contains an indefinite K = xx* − yy*, and four generic vectors give a rank-4 intensity map.
- (5) β ≤ 2 for m = (2,2,2): random three-dimensional subspaces violate the Knill–Laflamme condition. The general statement is the source's theorem.

Run: `python3 check.py` with numpy 2.4 and sympy 1.14 (exit 0; final line `bad= 0`; about 1 s).
