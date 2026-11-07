# Exact checks for Stein's finite-population sorting conjecture

Source: James D. Stein, *Properties of Subsets of Unobserved Items Constructed by Using Blackwell's Bet*, arXiv:2608.21955v1, Section IV's closing conjecture. The printed conjecture is E(U→A) ≤ E(U→B). Text extraction incorrectly drops the equality bar.

Run `python3 check.py`. The program uses only Python 3's standard library and exact `fractions.Fraction` arithmetic. Its output is in `check.out`; successful completion ends with `ALL_CHECKS_OK`.

A uniformly random labelled permutation gives disjoint examined/unexamined pairs. A uniform threshold assigns an unexamined value U to A with probability X, where X is its examined partner's value. The last item of an odd population is assigned fairly to A or B. Independent assignment randomness is integrated exactly. The expectations are ratios of expected assigned weight to expected assigned count, following the source's fractional-marble convention; they are not expectations of random within-run averages.

For N > 1, population mean m, population variance v, k = floor(N/2), and d = N mod 2, the program checks the pair identity E[XU] = m² − v/(N−1), the expected counts km+d/2 and k(1−m)+d/2, and the identities

- E_A = m − kv / ((N−1)(km+d/2));
- E_B = m + kv / ((N−1)(k(1−m)+d/2)).

It verifies all labelled permutations for 119 sorted populations on {1/4, 1/2, 3/4} of sizes 1 through 7 (204,555 permutations), plus 35 two-species populations with values 1/5 and 4/5 (46,232 permutations). Repeated values retain their labelled multiplicities. There are 70 odd grid populations and 21 constant grid populations. N = 1 is checked separately.

Fidelity checks reproduce Corollary 1.1's two-item formulas, giving E_A = 3/8 and E_B = 5/8 at S = 1/4 and L = 3/4, and the exact population mean 21/40 in the final simulation-table row (S = 2/25, L = 97/100, two copies each). The table's empirical conditional averages are not exact expectations and are not asserted as identities.

These are finite checks supporting a proposed general argument. They do not prove the unbounded conjecture and contain no Lean proof. The source's displayed Section IV recursion tracks assigned items without tracking removed examined items; these checks implement the disjoint-pair SAP in the prose, with the stated odd-singleton rule.
