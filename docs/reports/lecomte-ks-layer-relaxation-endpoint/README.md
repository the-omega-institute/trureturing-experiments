# Lecomte KS layer relaxation: attained endpoint

Source: Guillaume Lecomte, *The Certification Limits of KS-Type Layer
Relaxations: A Square-Root Ceiling and Its Breakdown*, arXiv:2607.12134v1,
Question 1, second sentence. Source license: CC BY 4.0. This directory contains
an independent check program and its output, with no copied source text.

Run with Python 3 and SymPy 1.14.0:

```sh
python3 check.py > check.out
```

The program uses integers, rational arithmetic and exact symbolic algebra.
It reproduces the source's three-atom energy and endpoint kernel constant,
and checks the chain-with-arm witnesses at populations 2 through 16.
It verifies the octahedron distance multiplicities and energy, the rational
zeta bounds, a coarse stability constant, and the predicted positive
endpoint-deficit bound 47/6606 > 1/141.

For the proposed family κ=0 and g(n)=γ−(C+K)(n−1), symbolic identities
check the validity comparison and strict exclusion threshold. Exhaustive
checks cover every ordered positive profile of size 1 through 12 at 81
parameter/reference-error settings (331,695 cases). A separate exact
equality case confirms that equality with the reference is not exclusion.
Fifteen actual oriented finite configurations check endpoint retention,
the singleton energy estimate and the discarded-to-retained interaction
estimate. The two-pile energy polynomial has quadratic coefficient
−κ(r) for any positive kernel entry.

These checks support the proposed classification F_LJ=(−C,−C1], with C1 the
supremum of binding energy per particle over oriented fully spanning
singleton-layer configurations. They do not determine C or C1, enumerate
all geometric configurations, prove validity for every configuration, or
prove the asymptotic classification. Those conclusions require the
general geometric and asymptotic arguments in the preregistration and
subsequent formal verification. This report contains no Lean proof.
