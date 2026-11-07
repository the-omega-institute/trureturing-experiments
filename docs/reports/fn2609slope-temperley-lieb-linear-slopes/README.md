# Spin-1 open-chain linear slopes

`check.py` builds the singlet projector and the antisymmetric triplet projector in the standard spin basis using exact rational entries. Run with Python 3 and SymPy:

```sh
python3 docs/reports/fn2609slope-temperley-lieb-linear-slopes/check.py
```

`check.out` contains the program's output. The program uses no floating point computations as evidence.

The checks reproduce the open-chain ground degeneracies and exact gaps for lengths 2–4, the length-3 slopes in Table 2 of Ferydouni–Nachtergaele, arXiv:2609.38074v1, and the length-4 gap rounded in that table. Ground-space triplet maxima are 1, 2, and (53+sqrt(457))/28 for lengths 2, 3, and 4. The program verifies their elementary bounds and the subadditivity inequalities whose total length is at most 4.

For lengths 3–5, it checks every insertion Gram block, the Hamiltonian insertion recurrence, the diagonal compression formula, and the local error bounds 2/3, 14/3 and 1. Rational Schur elimination certifies the operator norm inequalities on the compressed ground spaces. These imply the proposed error constant 36 after summation and normalization, subject to the general identities and estimates in the mathematical argument.

These finite checks support the proposed mechanism. They do not prove the limits as the chain length tends to infinity. No Lean proof is included.

Source: Robert Ferydouni and Bruno Nachtergaele, *On the Spectrum of Some Temperley-Lieb Spin Chains*, arXiv:2609.38074v1, Section 6.2.1 and Corollary 2.4, https://arxiv.org/abs/2609.38074v1.
