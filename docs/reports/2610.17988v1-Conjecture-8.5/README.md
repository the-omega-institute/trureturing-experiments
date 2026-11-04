# Clarke Conjecture 8.5: exact finite checks

Source: Oliver Clarke, *The equivariant Ehrhart theory of polytopes with
order-two symmetries: Part II*, arXiv:2610.17988v1, Sections 2 and 8.3.
The mathematical target is the unconditional coefficient comparison in
Conjecture 8.5. Its settlement is not proved by these experiments.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The stored output is
`check-output.txt`. A successful run exits 0 and ends with `ALL_CHECKS_OK`.

The program constructs supporting facets from rational nullspaces and uses
integer facet inequalities to enumerate lattice points. It uses no floating
point arithmetic or approximate hull routine. At dilation zero it uses the
source's singleton convention. An affine involution acts on a point of the
k-th dilate by `A*x + k*b`; its homogeneous matrix determines p and q.
The representation ring uses pairs of integers in the trivial/sign basis.
Orbit-parity assertions precede every multiplicity computation.

The fidelity checks reproduce the two reflected squares of Examples 2.3
and 2.4, the pyramid and bipyramid of Examples 5.1 and 5.2, and the
noneffective six-dimensional Example 5.6. In particular the last example
has h* = 1+t+t^3+t^4, H = 1+t^2+t^3+t^5 and sign multiplicity -1 in
degree 2. It satisfies Conjecture 8.5 while ordinary effectiveness fails.

The finite screen covers all 21 integer intervals with distinct endpoints
in [-3,3], all 102 full-dimensional invariant generating subsets of the
3-by-3 integer grid for reflection in x=0 or diagonal interchange, and
24 random full-dimensional three-dimensional polytopes (integer seed
261085). The sampled three-dimensional actions include affine reflections
with no fixed lattice point, simultaneous sign changes, and coordinate
interchanges. There are 141 cases with p>=q and six with p<q. All input
numerator coefficients through degree 2(d+1)+2 satisfy the comparison and
the normalization (1+t)^p phi(sigma) = (1+t)^q H.

For 56 construction steps the program separately enumerates the convex
hull B of P at height zero and the two exchanged apices (v,1) and
(sigma(v),-1). It checks the entire central slice, every noncentral slice,
fixed-point counts, multiplicities (p,q+1), and phi_B=(1+t)phi_P.
First steps use dilations k=0,...,4. Further balancing steps use k=0,1,2.
These are finite checks of the proposed geometric mechanism; the universal
slice bijection and the application of Clarke's Corollary 4.4 require a
mathematical proof and subsequent Lean development.

The experiment's enumeration is diagnostic only. The planned positive
settlement is a general convex-geometric and character-theoretic argument,
not a bounded positive enumeration. The general Ehrhart and equivariant
positivity infrastructure remains a formalization research line.
