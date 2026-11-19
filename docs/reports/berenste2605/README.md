# Commuting idempotents in type A 0-Hecke monoids

Run `python3 check.py` from this directory. The script uses integer permutations
and the right Demazure rule in §2.5 of Berenstein, Greenstein and Li,
*Hecke monoids, their homomorphisms and parabolicity*, arXiv:2605.09761v2.
It reproduces all three products and both alternation minima in Example 3.5.

The family checks use A_(r+2) with J = {1, …, r} and K = {r, r+2}
for 2 ≤ r ≤ 10. The script also checks every ordered support pair in A1
through A6 against the §3.1 expectation. The exact finite outputs are in
`check.out`.
