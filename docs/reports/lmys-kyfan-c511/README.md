# Lin–Miao–Yu–Sheng Conjecture 5.11: exact second-norm checks

Source: Z. Lin, L. Miao, G. Yu and H. Sheng, *On the Ky Fan k-norm of the
LI-matrix of graphs*, Electronic Journal of Linear Algebra 38 (2022), 680–696,
DOI: 10.13001/ela.2022.5943; arXiv:2012.11162v2.

Run `python3 docs/reports/lmys-kyfan-c511/check.py` from the repository root.
Requires Python 3 and SymPy (output supplied with SymPy 1.14.0). All matrix,
polynomial, ordering and norm calculations use exact integers and rationals.

The source-fidelity check reproduces the introduction's identity
`||LI(G)||_Fk = ||A(G)||_Fk` for regular graphs, for every k on K_n, 3 ≤ n ≤ 10.
The Figure-1 graph is rebuilt and its order, size, connectivity, pendant count
and Laplacian spectrum checked for all 28 admissible pairs with 3 ≤ n ≤ 9.

For 6 ≤ n ≤ 18, the separately constructed cones
`B_n = K_1 ∨ (P_3 ⊔ (n−4)K_1)` and
`H_n = K_1 ∨ (2K_2 ⊔ (n−5)K_1)` are connected, have n+1 edges, and satisfy
n ≤ n+1 ≤ 2n−3. Their different degree multisets rule out isomorphism.
Their characteristic polynomials give spectra
`{n,4,2,1^(n−4),0}` and `{n,3,3,1^(n−4),0}`. The program also checks the
characteristic polynomials of the centered matrices and their Gram matrices,
so the displayed singular-value and norm calculations have exact matrix checks.
Both second norms equal n. The third norms are different:
`n+2−2/n` and `n+1+2/n`.

These are finite experimental checks supporting an infinite-family refutation
outline. They do not constitute a Lean proof or an exhaustive extremal
classification. No third-norm uniqueness refutation is claimed. The graph uses
the source's intended n-vertex, m-edge construction, whose pendant count is
2n−m−3; it does not exploit the printed 2n−m+3 misprint.

Producer: codex-cli research flight 4d4a18fd9c0ea492bb09fcc5, thinking stage,
prior disclosure `repo-prior-exposed`. The source and literature screen are
reported in the separate preregistration draft, not asserted by this program.
