# Low-weight non-uniqueness converse

Source: M. Y. Avetisyan and R. L. Mkrtchyan, *Uniqueness of universal quantum
dimensions*, [arXiv:2608.30012v1](https://arxiv.org/abs/2608.30012v1), Section 3,
the unnumbered converse following Theorem 3.1. Equations (30)-(31) and Table 1
define the coefficient restriction and line weights.

Run with Python 3:

```sh
python3 check.py
```

`output.txt` contains the output. The program uses only integers and
`fractions.Fraction`; it performs no floating point calculations.

The fidelity checks reproduce the three line identities of equation (35)
and the four restricted multiset equalities asserted for equations (36)-(37).
The latter check tests exactly the multiset assertion, including any zero
arguments; it does not assert that the displayed source denominator restricts
to a nonzero function on those lines.

For every one of the 160 subsets of total weight below three, the construction
uses even and odd subset sums shifted by (1000, 1000000, 0). The empty selection
uses the auxiliary vector (1, 0, 0). The check verifies the coefficient bounds,
the subset-toggle pairing and restricted multisets, and a real point on each
of the twelve lines where every numerator and denominator argument is nonzero.
At x=4 and (alpha,beta,gamma)=(1000000,-1000,1000000000), the empty subset gives
a zero numerator argument and all denominator arguments are positive. The
minimum denominator argument is 998000, and there are at most 16 factors on
each side.

The transcendental implications use the standard real facts that sinh(t)=0
exactly when t=0, and sinh(t)>0 when t>0. The script certifies their exact
arithmetic premises and coefficient identities. It is corroboration for a
uniform subset-toggle proof, not a Lean proof or a numerical approximation
to one. No proof of the converse or its negation in Lean is included.
