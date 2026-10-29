# Simple evolution algebra automorphism maximum

Exact finite checks for Conjecture 4 of Ivan Kaygorodov, Cándido Martín
González and Antonio Viruel, *Evolution algebras: a survey of results and a
list of open problems*, arXiv:2610.22248v1, Section 2.3.

Run `python3 docs/reports/kmv-c4-simpleaut/check.py` with Python 3. The script
uses only the standard library: integer arithmetic and `fractions.Fraction`.
The checked output is `check.out`; its final line is `ALL_CHECKS_OK`.

- **Source fidelity:** enumerate every labeled support digraph, allowing
  loops, through dimension four. Keep exactly the strongly connected graphs
  admitting a support permutation, a necessary condition for an invertible
  structure matrix. The respective support counts are 1, 4, 138 and 24,416.
  The maxima of diagonal order times full graph group order are 1, 6, 21
  and 60, reproducing the source's stated verification through dimension four.
- **Independent diagonal check:** eliminate scalar assignments along actual
  directed paths and take a gcd of the remaining integer exponents. For all
  143 retained supports through dimension three, independently compute the
  order from the gcd of all maximal minors of the integer relation matrix.
- **Proposed mechanism:** check the directed-cycle period, every edge's
  increment, equal sizes of its classes, the rotation action of every graph
  automorphism, and the joint diagonal/permutation bound on every retained
  support through dimension four.
- **Attainment:** for dimensions one through eight, compute exact
  determinants for the cyclic matrix and `I + J`, enumerate all graph and
  coefficient-preserving permutations, and check every cyclic root assignment
  as integer phases. Together with monomiality for perfect evolution algebras,
  these give full group orders `n(2^n - 1)` and `n!` in that range.
- **Arithmetic:** check all 2,191 pairs `h,k >= 2` with `hk <= 500` in the
  factorial comparison, and all 4,096 pairs `1 <= a,b <= 64` in the Mersenne
  gcd identity.
- **Controls:** detect a strongly connected singular support with unequal
  period classes and an invertible weighted algebra whose graph transposition
  does not lift. No surjectivity onto the full graph group is assumed.

These are finite checks, not a proof for unbounded dimension. The graph
calculation bounds the actual permutation image by the full graph group;
it does not assert that every symmetry lifts for arbitrary weights.
