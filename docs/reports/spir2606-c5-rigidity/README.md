# Exact checks for C5 resonance rigidity

Source: Călin Spiridon, *On monomial resonance*, arXiv:2606.19006v1,
Question 3.29 and Proposition 3.28.

Run with Python 3.12 or later and SymPy 1.14:

```sh
python3 check.py
```

`output.txt` contains the complete output. The program uses integers, rational
row reduction with `fractions.Fraction`, and exact symbolic polynomial algebra.
It exits zero and prints `ALL_CHECKS_OK` only after every assertion passes.

The checks reproduce the graph-resonance formula on every nonempty support of
all 64 labelled graphs on four vertices and all 1024 labelled graphs on five
vertices (32,704 supports). They recover the two components for C4 and the five
components for C5. On coordinate subspaces, only the original cycle has the
same resonance as each respective cycle.

For source fidelity, the program also checks Proposition 3.28 on 140 seeded
random rational subspace presentations and 20 changed bases of the C4
subspace. Every unequal sampled subspace has an explicit rational vector with
a different contraction-rank resonance test. The twenty equality cases are
certified by equality of rational row spaces, rather than grid agreement.
The decomposable-bivector identity in the source proof is checked symbolically.

For the proposed C5 argument, it checks pair incidence, the five nonzero integer
determinants showing a common resonance-free two-plane, the universal pencil
factorization and homogeneous cubic determinant, and all 203 possible patterns
of zero component restrictions and coincidences of the remaining lines.

These checks support a structural paper argument about arbitrary subspaces
over algebraically closed characteristic-zero fields. They are neither an
exhaustive search over arbitrary subspaces nor a Lean proof. The root-existence
step and the passage from resonance components to restriction ranks are
mathematical arguments outside the program.

The check program is original work; its SymPy dependency is distributed under
the BSD-3-Clause license. No third-party program source is copied here.
