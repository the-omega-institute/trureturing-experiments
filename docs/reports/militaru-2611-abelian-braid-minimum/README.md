# Minimum braid carriers for finite abelian groups

`check.py` uses Python 3 and the standard library only; all evidence uses integers,
finite sets, tuples, and exact permutation operations. Run:

```sh
python3 docs/reports/militaru-2611-abelian-braid-minimum/check.py
```

`output.txt` is the program's output. It ends with `ALL_CHECKS_OK`.

The fidelity check reproduces Militaru, *Centralizers of maps as automorphism
groups of Yang–Baxter solutions*, arXiv:2611.03380v1, Example 3.3 and Remark 3.4:
`R(x,y)=(y,x+1)` on the cyclic group of order n has full automorphism group Cn.
It checks n from 2 through 8, enumerating every carrier permutation.

The proposed construction has blocks of primary cyclic orders q_i. Its left
translation on block i translates every block j with i ≤ j by one; the associated
braid map is `B(x,y)=(lambda_x(y),x)`. The program checks the braid equation on
every triple, bijectivity, and both nondegeneracy conditions. It enumerates all
carrier permutations and proves equality of the resulting finite set with the
explicit block translations. It checks that the parameterization is bijective,
recovers its inverse, and preserves the complete group multiplication table.
Thus the finite isomorphism check is stronger than comparing element orders.

The eight cases are C2², C2×C4, C4×C2, C3², C6, C2³, C4² and C2²×C3.
For C6 the primary factors are 2 and 3, so the sharp carrier has five points;
the cyclic six-point fidelity example is a separate, nonoptimal construction.
The Chinese remainder isomorphism is checked explicitly.

For each group the program enumerates every subgroup, computes the intrinsic
primary-factor weight from the p-torsion filtration, checks subgroup monotonicity,
product additivity and the quotient weight bound, and solves the minimum-cost
kernel-intersection problem over all coset-orbit families. Its minimum faithful
permutation degrees are respectively 4, 6, 6, 6, 5, 6, 8 and 7.

These are finite consistency checks. The proposed universal result requires the
general triangular automorphism argument and the classical finite abelian minimal
faithful-degree theorem. The latter is stated on p. 207 of Easdown–Praeger,
*On minimal faithful permutation representations of finite groups*, Bull. Austral.
Math. Soc. 38 (1988), 207–220, DOI 10.1017/S0004972700027489, with attribution
to Johnson (1971) and Karpilovsky (1970).
