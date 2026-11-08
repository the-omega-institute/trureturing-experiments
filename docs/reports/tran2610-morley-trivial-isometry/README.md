# A simplex with regular Morley simplex and trivial full isometry group

The check certifies a positive answer to the first clause of Quang Hung Tran,
*Morley Simplices in Higher Dimensions: Regularity, Reflections, and
Counterexamples*, arXiv:2610.01216v1, section 12, item 4, in dimension 2499.
It also reproduces the source's certified Johnson example for k = 9
(Theorem 10.2): spectral compatibility changes sign on the stated interval
p in [1.6676, 1.6680], with the required positive normal eigenvalues.

Run with Python 3 and SymPy (BSD licence):

```sh
python3 check.py
```

The stored run uses Python 3.14 and SymPy 1.14.0. All evidence uses integers,
exact Fraction interval endpoints, or exact symbolic algebra. Printed decimals
are terminating rationals rounded outwards to 36 decimal places. The final
line is `ALL_CHECKS_OK`.

The explicit order-50 Latin square is included in check.py. Its graph is checked
exhaustively to have parameters (2500,147,50,6). Every collinear triple has
exactly its line as its common closed neighbourhood. A five-clique contains
a collinear triple, by the pigeonhole principle applied to the three line
types through one vertex. Hence the maximum 50-cliques are exactly the 150
rows, columns and symbol lines. Their disjointness identifies the three
parallel classes intrinsically. The cycle partitions between disjoint lines
provide 150 pairwise distinct line signatures. Every possible line image
under a graph automorphism is checked and forced to be itself. Intersections
of lines then force every graph vertex to be fixed. This certifies the full
graph automorphism group, including permutations of the three coordinate
classes, without an external graph-automorphism program.

For epsilon = 1/50 and phi(x) = (3x-x^3)/2, set a = phi(p) and
b = -(epsilon^2 + 3 epsilon (1-epsilon) a)/((1-epsilon)(1-2epsilon)).
There is a unique r in (1,2) with phi(r) = b. The continuous compatibility
function F has opposite strictly certified signs at p = 8/5 and p = 13/8.
Rational bisection retains an interval containing an actual zero by the
intermediate value theorem; it does not assume root uniqueness. The check
certifies 1 < p < r < 2, positive normal eigenvalues and negative nonconstant
K eigenvalues throughout its final enclosure.

Writing B = 1+2b-3a, C = 1+2r-3p, H = I+aA+b(J-I-A),
K = I+pA+r(J-I-A), and P = I-J/2500, exact symbolic algebra verifies
C^2 H - B K P K = F(-150 I-50 A+3 J). At an actual zero this gives
H = tau K P K with tau = B/C^2 > 0 and H*1 = 0.
On the Euclidean space V = 1-perp, set u_i = sqrt(tau) K P e_i and
Delta = {x: 1+<u_i,x> >= 0 for all i}. The affine map
x -> (1+sqrt(tau)Kx)/2500 identifies Delta with the standard simplex.
With t = 2500/kappa_0 the points M_i = (t/sqrt(tau)) P e_i satisfy the
source's nearer-trisector equations and all have pairwise squared distances
2t^2/tau. Any self-isometry fixes the unique incenter and preserves H,
so it preserves A; conversely graph automorphisms act isometrically.
Thus Delta has trivial full Euclidean isometry group.

The finite graph certificate supplies an input to this general geometric
argument. The record contains no Lean proof and supplies no kernel-verified
settlement. The claim concerning distinct edge lengths, the dimension-4/5
reflection conjectures, and existence in every dimension are outside this
result. An asymptotic Latin-square existence theorem is not used.
