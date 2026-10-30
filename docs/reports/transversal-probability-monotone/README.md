# Common transversal probability under normal extension

`check.py` checks Aivazidis–Loukaki–Müller, *On the common transversal
probability in finite groups*, arXiv:2209.00986v2, Problem 5.2.

Run with Python 3 and its standard library:

```sh
python3 docs/reports/transversal-probability-monotone/check.py
```

For every subgroup H and every normal subgroup N of each constructed group,
it checks P_G(H) ≤ P_G(HN), and equality precisely when N ≤ H or H is normal.
It checks the quotient-image identity and the orbit/fiber comparison on actual
cosets. Probabilities are obtained independently from the weighted permanent
of the actual left/right coset intersection matrix, then compared with the
factorial product in Proposition 3.2. All arithmetic uses integers and
`fractions.Fraction`.

The group list contains the standard 24 isomorphism types of order at most
12: the cyclic groups of orders 1–12, C2×C2, S3, C4×C2, C2×C2×C2, D4, Q8,
C3×C3, D5, C6×C2, D6, A4 and Dic12. Here Dn has order 2n, and Dic12 is given
by a^6=1, b^2=a^3, bab⁻¹=a⁻¹. Additional groups are D7, D8, S4 and S3×S3.
Each multiplication table is checked for identity, inverses and associativity;
every subgroup is enumerated by adjoining generators until closure.

Fidelity checks reproduce Lemma 4.2(ii) for D3–D8 and Lemma 4.2(iv),
tp(A4)=2/9, as well as Corollary 3.4. The stored output reports 1,807 pairs,
693 strict inequalities, and ends in `ALL_CHECKS_OK`.

This is finite computational evidence. It does not prove the universal
statement or certify that the problem remains open in the literature.
