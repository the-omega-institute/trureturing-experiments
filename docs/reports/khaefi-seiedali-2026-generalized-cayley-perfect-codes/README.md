# Generalized Cayley subgroup regular-set check

`check.py` checks the proposed extension of Theorem 4.7 from arXiv:2611.00893v1. The input `smallgroups.dat` was emitted by GAP 4.16.1 with the SmallGrp package using `SmallGroup(n, i)` for every group order `n ≤ 16`, every involutory automorphism, and every proper invariant subgroup. The Python check uses integer multiplication tables, bit masks, and exact finite-set cardinalities. It tests all admissible internal and external degrees, the subgroup-perfect-code obstruction, the cyclic-index-two source family, the trivial subgroup cases, and involution-fixed right cosets.

Commands:

```text
/Users/auric/.op-scratch/op/ks2611p2/gap/gap -q < gap-export.g > smallgroups.dat
python3 check.py > check-output.txt
```

The final output line is `ALL_CHECKS_OK`.
