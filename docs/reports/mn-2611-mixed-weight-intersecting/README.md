# Mixed multiplicative weights for maximal intersecting families

The source is Elham Mehdi-Nezhad, *Domination and maximal intersecting families in comaximal graphs of commutative rings*, arXiv:2611.08022v1, Problem 7.2. This report tests a proposed answer; it is not a formal proof or a novelty certificate.

Run `python3 check.py`. The program uses only the Python standard library, and all numerical comparisons use exact integers and `fractions.Fraction`. The stored output is `output.txt`.

The independent enumerator is Bron–Kerbosch on the intersection graph of all nonempty proper subsets. It finds 2, 4, 12, 81 and 2646 maximal intersecting families for ground-set sizes 2 through 6, and checks their defining intersection and maximality conditions directly.

For positive weights `w(A) = c product(t_j, j in A)`, let `L = {j : t_j < 1}`, let `H` be its complement, and write `t_B` for the product over `B`. When `L` is nonempty the proposed minimum is

    (c/2) product(1+t_j, j in H) sum(min(t_B,t_(L\B)), B subset L) - c t_U.

The checks compare this expression with the actual minimum over every enumerated family, verify every proposed tie-breaking minimizer, and check independence from strictly expensive coordinates for every optimizer in these finite instances. They cover 132 independently seeded mixed vectors for sizes 2 through 5, another 33 for size 6, and the complete 117-vector grid over `(1/2,1,2)` for sizes 2 through 4. The orchestrator's original random vectors were not supplied; this is an independently seeded reproduction of its counts and test size.

Fidelity checks reproduce Remark 4.5 (three family weights 10 and one 12, against complementary-pair lower bound 6), the Boolean-weight formula, and the proposed nondegenerate mixed example with minimum 19. Eighty further unmixed vectors check agreement with Theorem 4.4 after its stated reparametrization and with Proposition 6.5. The mixed formula overlaps Theorem 4.4 whenever all parameters are at most 1 and at least one is strictly smaller. It does not apply to strictly greater-than-1 vectors. Its extension to empty `L` agrees with Proposition 6.5 on the boundary where the smallest parameter equals 1; the proposition supplies the value for the remaining all-high vectors. No mixed-domain assertion is inferred from those disjoint-domain comparisons.

Finite checks do not establish the formula for arbitrary real parameters or unbounded ground-set size. Those require the proposed general argument using increasing self-dual cube families and coordinate bias normalization.
