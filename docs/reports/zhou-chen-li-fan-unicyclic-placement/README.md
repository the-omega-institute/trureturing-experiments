# Zhou–Chen–Li–Fan unicyclic placement

Exact independent checks for Conjecture 1, odd branch, in Yu Zhou, Yuanyuan Chen, Dan Li and Dandan Fan, *Minimizing the spectral radius of unicyclic graphs with fixed order and dissociation number*, arXiv:2610.38539v1, Section 5.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The saved output is `check.out`. All evidence uses integer or rational arithmetic and exact symbolic algebra.

The fidelity check reproduces the four-vertex triangle-plus-leaf example in the proof of Theorem 1: its dissociation number is 3. Exhaustive binary cycle checks also give dissociation number floor(2m/3) for 3 ≤ m ≤ 18.

For 0 ≤ a ≤ 6, the two explicitly constructed unicyclic graphs have n = 12a + 23 and dissociation number phi = 12a + 17. The upper bound restricts any dissociation set to the induced C18; deleting the six hubs attains it. All allowed low-hub leaf placements of the opposite-high graph are equivalent under cycle automorphisms.

Symbolic checks verify the two-subdivision and pendant elimination equations, the reduced determinant difference, and the reflection identity supporting the proposed infinite-family argument. These identities and finite checks are not a formal proof of the universal spectral comparison.

Independently, exact positive integer vectors on the actual 23-vertex graphs certify rho(B_0) < 213/100 < rho(A_0). Every edge, vector coordinate and integer residual appears in the output. The final line is ALL_CHECKS_OK.

The intended source reading restores the missing minus one in the odd-case high-multiplicity numerator and interprets the balancing index as the next selected hub. The preregistration quotes the printed TeX without these corrections.
