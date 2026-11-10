# Fixed pins on lacunary integer spheres

`check.py` checks finite instances of the family proposed to refute Palsson–Wang, *On Densities in the Erdős-Falconer Distance Problems over the Integers*, arXiv:2610.10665v2, Conjecture 8. It uses Python's standard library and exact integers and rational numbers throughout. Run `python3 check.py`; the expected final line is `ALL_CHECKS_OK`. The captured result is `output.txt`.

Fidelity checks reproduce Example 5.1: every integer pin of the one-dimensional lattice gives the perfect squares, and the four-dimensional lattice realizes all squared distances through the four-square theorem. These are finite checks of the explicit descriptions, not numerical estimates of their dimensions.

For each tested dimension, the program builds projected boxes on spheres using the lexicographically first nonnegative four-square completion. It verifies shell sizes, projection injectivity, exact norms, the finite affine basis, rational shell masses, geometric partial sums, sampled cube mass estimates, Cauchy–Schwarz and distance identities, shell distance-count bounds, the lower count at the pin e₁, and the strict threshold for dimensions at least nine. Shells are tested for dimensions 5 through 12, with shell indices at most three and at most 40,000 points per shell. Every point in each selected shell is checked; only the cube and pin sets are sampled.

The proposed infinite construction has discrete Hausdorff dimension d−4 and maximal pinned upper discrete mass dimension 1/2 for d≥9. Finite computation does not establish those infinite statements. Their covering and counting arguments belong to the preregistration and subsequent formal proof. Dimensions 5 through 8 do not yield a strict contradiction from this family's d−4 exponent.

Source definitions: arXiv:2610.10665v2, Definitions 2.1, 2.3 and 2.4 and Example 5.1. The lower cutoff for lattice cubes is Bishop–Peres, *Fractals in Probability and Analysis*, Definition 8.1.2, printed p.245. Existing mathematical ingredient: Lagrange's four-square theorem; Mathlib supplier `Nat.sum_four_squares`.
