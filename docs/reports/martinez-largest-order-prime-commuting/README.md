# Martínez largest-order-prime commuting bound

Source: Juan Martínez, *On the commuting probability of π-elements in finite groups*, [arXiv:2307.09868v3](https://arxiv.org/abs/2307.09868v3), Theorem A and the unnumbered largest-prime question at the end of §4; DOI [10.1002/mana.202300338](https://doi.org/10.1002/mana.202300338).

Run with Python 3 and the standard library:

```sh
python3 check.py
```

The output is in `check.out`. Every count uses integers and every ratio uses `fractions.Fraction`.

The fidelity check constructs PSL(2,5) directly as determinant-one matrices modulo 5, identified modulo ±I. It reproduces the source's formula Pr₅(PSL(2,5)) = 29/125, from 25 elements and 145 commuting ordered pairs. It also checks equality in Theorem A for a transposition in S₃ with π = {2}.

For each prime q in {5,7,11,13,17,19,23,29,31}, the program constructs G = S₃ × D₂q and x = ((12),r), with π = {2,q}. It counts |G| = 12q, o(x) = 2q, |Gπ| = 8q and |CG(x)π| = 2q, checks the trivial center, and computes Oπ(G) = 1 × D₂q. Hence the ratio is 1/4 > 1/q and x is outside the radical. It also checks the source's 1/2 bound for every eligible element in these groups.

The tables of each factor and of PSL(2,5) are checked for associativity, identity and inverses. Direct products use coordinatewise multiplication. Radicals are computed independently of the predicted product formula: x belongs to a normal π-subgroup exactly when its normal closure is a π-group. The program generates those normal closures from all conjugates, and verifies that the resulting radical is a normal π-subgroup. This certifies maximality because any normal subgroup containing x contains its normal closure.

These finite computations support the proposed counterexample family. The argument for every prime q ≥ 5 and a Lean proof remain separate obligations. The underlying products of odd-prime dihedral groups already occur in Detomi–Senise–Shumyatsky, [arXiv:2605.22955v1](https://arxiv.org/abs/2605.22955v1), Example 3.3, for Sylow-average and Fitting-index questions; no novelty is claimed for that construction.
