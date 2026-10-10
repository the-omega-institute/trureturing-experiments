# Cruz Ortega–Menéndez-Conde: the resistance-constrained Fiedler maximum is not attained at invariant weights

`check.py` supports a refutation of the conjecture of arXiv:2306.08187v1 (abstract and Section 5): "for every graph the maximum is attained for weights that are invariant under automorphisms". The maximum in question is that of λ₂(G; w) subject to R(G; w) = R(G), where R is the total effective resistance over all unordered vertex pairs.

The counterexample family is K₂,ₘ for m ≥ 3, with hubs u, v and leaves x₁, …, xₘ.
- With conductance 2 on each u–xᵢ edge and 1 on each v–xᵢ edge (Laplacian Bₘ), λ₂ = 3 and R = (2m³+2m²+5m+9)/(6m).
- The unweighted graph has λ₂ = 2 and R = Dₘ/(2m), where Dₘ = m³+m²+2m+4.
- Scaling Bₘ by cₘ = (2m³+2m²+5m+9)/(3Dₘ) meets the resistance budget, with λ₂ = 2 + (m+1)/Dₘ > 2.
- The automorphism group acts transitively on the edges, so every invariant feasible assignment is the unweighted one, which has λ₂ = 2.

What is checked, exactly with sympy, for m = 3..8:
- the Laplacian spectra;
- the Kirchhoff index through the pseudoinverse;
- the source's harmonic-current resistance for each pair class;
- the scaled assignment.

Edge-transitivity of the automorphism group is checked by enumeration for m = 3, 4. For m = 3 the scaled assignment has λ₂ = 48/23 and R = 23/3.

Run: `python3 check.py` with sympy 1.14 (exit 0; final line `bad= 0`).
