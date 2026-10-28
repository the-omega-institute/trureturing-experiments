# Quantum Černý extremal-word checks

Source: Pui Hang Lee, Pui Yee Lee and Bjørn Kjos-Hanssen, *Quantum Černý complexity of binary words*, arXiv:2609.40154v1, Theorem 3.1 and §6 open problem 4, https://arxiv.org/abs/2609.40154v1.

Run with Python 3 and SymPy 1.14.0:

```sh
python3 check.py > output.txt
```

The script uses integers and exact SymPy algebra, including rational numbers and square roots. It uses no floating point or numerical eigenvalue thresholds.

For fidelity it reproduces the source's constant-word nilpotent model for lengths 1–12: the ranks decrease by one, the only zero word of length at most m is 0^m, and the dimension threshold is ceil(sqrt(m+1)). For m=8 it constructs the corresponding d=3 quantum channels and verifies the actual density matrix and Choi matrices exactly, reproducing the construction underlying qc(0^8)=3.

For the proposed refutation family 0^(m−1)1 it enumerates every binary word of length at most m for 2≤m≤12. Only the nominated word annihilates the entire m-dimensional model; all other such words keep its initial basis vector nonzero. It checks the descending image ranks m,m−1,…,1,1. There are 32,733 word tests across both families.

For 2≤m≤8 it realizes the marker model as actual d-by-d channels with d=ceil(sqrt(m+1)), ε=1/[2d(m+1)], using orthonormal traceless Hermitian matrices. It checks trace preservation, unitality, the images of every selected basis matrix, and the reachable states that span the model. Exact LDL* factorizations have strictly positive pivots and reconstruct the density matrix and J(Φ)−I/(2d), with the unnormalized Choi convention. This certifies complete positivity with the stated margin. Frobenius-norm identities are checked separately.

These are finite checks of the proposed argument. They do not enumerate all quantum channels or prove the universal lower bound, all-length equality, literature novelty, or the eventual-uniqueness refutation in Lean. The universal image-chain argument and the paper's general upper-bound construction belong to the preregistration's proof outline.

The expected final output line is `ALL_CHECKS_OK`.
