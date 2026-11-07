# Equational spectrum and abelianness

Exact preregistration checks for Carles Cardó, Keith A. Kearnes and Ágnes
Szendrei, *Finite algebras whose equational probability spectrum contains a big
gap*, [arXiv:2610.37841v1](https://arxiv.org/abs/2610.37841v1), Section 4,
question (4)(a). This report supplies computational evidence for a proposed
negative answer; it is not a Lean proof.

Run with Python 3, using only the standard library:

```sh
python3 docs/reports/cardo-kearnes-szendrei-spectrum-abelianness/check.py
```

All probabilities, counts, slopes and inverse interpolation checks use integers
or `fractions.Fraction`. `output.txt` is the script's complete output and ends
in `ALL_CHECKS_OK`.

The source fidelity check reproduces Example 4.1's semilattice equation formula
for 57 triples with a,b,c in 0 through 3, subject to a+c>0 and b+c>0.

For the semigroup S with rows and columns ordered 0,a,b and multiplication table
[[0,0,0],[0,b,0],[0,0,0]], the script constructs term functions by closing the
projections under the actual binary operation in arities 1 through 6. The result
matches the projections, squares, products of two distinct variables, and zero.
Every comparison probability is in

Q = {1/3,5/9,17/27,2/3,7/9,65/81,23/27,8/9,1},

and all nine values occur from arity 4 onward. It checks all eleven comparison
types and the term-condition failure t(0,0)=t(0,a), t(a,0) != t(a,a).

For cyclic additive groups the script checks exact uniform fibers for all unary
coefficient residues modulo 3^e for e<=5, selected residues for 6<=e<=8, and all
binary coefficient vectors modulo 3,9,27. All nine values 3^(-j), 0<=j<=8, are
attained by native equations (3^(8-j)+1)x=x modulo 6561. It also checks 532170
linear term-condition instances modulo 3 and 9 and 190 direct-product equation
counts on S^2 and S^3.

The proposed family is A_d=S^(8d), B_d=(Z/6561Z,+)^d for every d>=1. These
algebras share cardinality 3^(8d) and the binary signature. The proposed spectra
are {q^(8d):q in Q} and {3^(-jd):0<=j<=8}. The script checks d=1 through 6:
nine distinct ordered spectral points, the common minimum 3^(-8d), positive
rational slopes of the interpolation fixing 0 and 1, its inverse at interval
midpoints, and the lifted nonabelian witness. At d=1 both algebras have size 6561.

The general argument still requires arbitrary-arity normal-form induction,
uniform fibers of homogeneous linear congruences, the native power-spectrum
law, additive cancellation for the term condition, and the piecewise linear
homeomorphism construction for all d. Finite checks alone do not establish these
universal statements. Nonabelian finite-spectrum semigroups are already covered
by the source's Example 4.8; the proposed contribution is the exact matching
family. The power-spectrum law is already in Cardó,
[arXiv:2604.08446v1](https://arxiv.org/html/2604.08446v1), Section 3.
