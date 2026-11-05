# Multiorbital versus orbital bivariate chromatic polynomials

Source: Melanie Gerling, *The Multiorbital Bivariate Chromatic Polynomial:
A Subgroup-Lattice Refinement*, arXiv:2609.01866v4, 15 October 2026,
Introduction/Conclusion converse question, Definition 3.8 and Proposition 3.9.

Run `python3 check.py` with Python 3 from this directory. Only the standard
library is required. `check.out` is its complete output; success is exit 0
and the final line `ALL_CHECKS_OK`. Every numerical or polynomial comparison
uses integers or `fractions.Fraction`.

The source fidelity check reconstructs the Section 3.3 examples: the reflected
three-vertex path and an edge plus an isolated vertex with the edge swapped.
It retains the loop in the second quotient, obtains their common orbital
polynomial `(x^3+x^2)/2-xy`, and checks that their multiorbital difference is
`y(1-x)`. Inclusion-exclusion and direct coloring enumeration are compared
on every quotient for `0 <= y <= x <= 4`. The star formula
`S_m=(x-y)x^m+y(x-1)^m` is independently checked by inclusion-exclusion and
literal coloring counts for `m=0,...,5` and the same palettes.

For the proposed converse refutation use `G=C2 x C2 x C3`, with elements
`(u,v,t)`. The transitive blocks are `T=C3` with translation by `t`, `A=C2`
with translation by `v`, `B=C2` with translation by `u`, `C=C2` with translation
by `u+v`, the regular block `R=G`, and a fixed singleton `I`. The two leaf sets
are `4T+9A+9B+5C` and `R+8A+8B+14I`, both of size 58. Add a fixed center
adjacent to every leaf. These are connected simple stars with faithful actions.

The check enumerates all 4096 subsets of G and finds exactly ten subgroups,
checks their product decomposition, and computes each coefficient
`c_G(g)=sum(1/|H| for H containing g)`. It constructs the actual permutations,
checks every one of the 144 composition pairs for each action, and checks
faithfulness, graph preservation, all twelve quotient graphs and cycle counts.
It then compares the subgroup-average and element-weighted formulas as sparse
polynomials in two variables with rational coefficients, giving

```text
12F1=12F2=44S58+22S50+24S44+12S40+12S36+6S32
OP2-OP1=(S40-S36)/4
OP1(2,1)=24209895276412929
OP2(2,1)=24210152974450689
OP2(2,1)-OP1(2,1)=257698037760
```

Adding r fixed leaves increases every leaf-cycle count by r, for every group
element. The checked equality of weighted exponent profiles therefore persists
under the symbolic replacement `m -> m+r`. The difference becomes
`(S_(40+r)-S_(36+r))/4`, whose value at `(2,1)` is `15*2^(34+r)>0` for every
integer `r>=0`. The program additionally constructs and checks the actions and
expanded polynomials for `r=1,2,10,100`. These additional samples are finite;
the statement for all r uses the stated singleton-cycle shift argument.

This report is exact computational evidence and a mathematical proposal.
It contains no Lean proof and does not claim a kernel-checked settlement.
