# Nearly Gorenstein G-polytopes and levelness

`check.py` checks Hibi--Seyed Fakhari, *Bounded powers of edge ideals: Nearly
Gorenstein polytopes*, arXiv:2610.24438v1, Conjecture 0.2, using integers only.
Run with Python 3.12 or newer:

```sh
python3 docs/reports/hsf2610-nearly-gorenstein-level/check.py
```

The output enumerates 121 cases: every labelled graph without isolated vertices
on three vertices with all capacities in {1,2}; every such graph on four
vertices with capacities (1,1,1,1) and (1,2,1,2); and seven additional cases.
The latter include C5, unequal capacities on C5 and a triangle with a pendant
edge, a disconnected graph, codegree-one examples, and source Example 2.10(3).
There are 111 nearly Gorenstein cases, including 65 that are not Gorenstein.
All 111 are level. One of the remaining cases is not level.

## Construction and exact decisions

Repeated products of edge exponent vectors, truncated by the capacities, give
the last nonempty bounded power and its minimal generators B. Taking all
coordinatewise divisors gives D. The program checks downward closure and the
augmentation axiom directly, and checks that D contains 0 and every unit vector.
It computes rho(A)=max{x(A):x in D}, verifies monotonicity and submodularity,
and independently checks that D is exactly the nonnegative integer solutions
to x(A)<=rho(A). The classical integral-polymatroid representation and normality
theorems identify this polytope with conv(D) and its semigroup ring with the
lattice points of its homogenized cone.

Degree-j canonical monomials are enumerated by x_i>0 and x(A)<j rho(A) for
every nonempty A. A monomial is a minimal canonical generator precisely when
no subtraction of a degree-one lattice point leaves a canonical monomial.
All canonical generators have degree at most n+1. To see this bound, express
an interior point of the (n+1)-dimensional cone as a nonnegative combination
of at most n+1 height-one vertex rays. Its degree is the sum of the coefficients.
Above degree n+1, a coefficient exceeds one; subtract that ray. Every previously
positive coefficient stays positive, so every strict cone-facet inequality
stays strict. This gives a lower-degree interior lattice point. Induction
proves the asserted generator bound. This argument needs no levelness assumption.

Thus the finite generator list determines the entire canonical ideal. An
inverse monomial v has v+g in the ring for every canonical generator g.
Testing this condition exactly decides inverse membership. A degree-one
monomial a lies in the trace iff a=g+v for some canonical generator g and
inverse monomial v. Because the ring is generated in degree one, this decides
the source's whole irrelevant-ideal containment. Levelness is independently
decided by whether all minimal canonical generators have the codegree.

For every nearly Gorenstein case, the program also checks the proposed exact
lattice sum U_N=U_s+((N-s)P intersect Z^n) through N=n+1. For the non-Gorenstein
cases, every tested inverse factorization of a degree-one monomial uses
canonical degree s. This checks the degree restriction needed in the proof.

## Source fidelity

For all tested complete graphs, the program reproduces the source's complete
classification in Theorems 1.5 and 1.6 and hence Corollary 1.7: nearly Gorenstein
implies level. It also reproduces Example 2.10(3): K_(2,2), capacities (3,2,2,2),
is nearly Gorenstein with codegree one (and has two canonical generators).

These are finite exact checks with the stated classical semigroup and convex
geometry identifications. They provide no Lean certificate and no proof of
the conjecture for arbitrary graphs or capacities. The proposed positive
settlement uses a general mathematical argument, not finite enumeration.
