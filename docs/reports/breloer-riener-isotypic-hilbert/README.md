# Breloer–Riener isotypic Hilbert series: exact checks

Source: Henri Breloer and Cordian Riener, *Symmetric and Isotypic Hilbert Series for Symmetric Ideals*, arXiv:2606.18895v2, Example 5.11 and Conjecture 5.12.

Run `python3 check.py`. It uses only Python's standard library and integer arithmetic. The stored output ends with `ALL_CHECKS_OK`.

The fidelity check independently enumerates multisets of local exponent vectors to reproduce all three series in Example 5.11 through degree 9. It also checks the third example's stabilization at two sites.

For each r = 2, 3, 4, 5, the family check enumerates the normal monomials of the local quotient with variables z, x₁, …, xᵣ and relations xₐxᵦ = 0 for a < b. It verifies local degree counts rm + 1 and stable orbit counts through degree 7, and HᵣP₁^(r−1) = P₂^r through degree 16. It checks the plane equations on a finite integer grid and the exponent separation used after clearing P₁ denominators.

These are finite checks. The general local basis, orbit correspondence, geometric component restriction, analytic convergence, and algebraic independence arguments belong to the preregistration's mathematical outline. The program does not certify infinite-series nonmembership or prove the conjecture's negation.
