# Zhao fixed-dimension Recht–Ré optimum: exact checks

Source: Jinze Zhao, *Sharpness and Rigidity of the Two-Factor Recht–Ré Lower Bound*, arXiv:2610.00155v1, Sections 2–4, Corollary 6.1 and the first remaining question in Section 8.

Run with Python 3 and Sympy 1.14.0:

```sh
python3 check.py
```

`output.txt` is the deterministic output. All matrix entries and comparisons use integers, rational numbers, or exact algebraic numbers. No floating-point computation is used.

The checker reproduces the source's canonical extremizer identities and normalization for n=2,…,9. It checks the proposed attained value −(k−1)/(4k), k=min(n,d), on 64 pairs with 2≤n≤9 and 1≤d≤8, including exact unit eigenvectors and positive projection certificates. It also checks 128 positive definite perturbations.

For 60 reproducibly generated rational feasible tuples it verifies every PSD principal-minor condition, the ordered-sum identity, 55 positive-pair merges, 59 removed annihilators, the resulting quadratic-form comparisons and the source's slack identity. These tuples are samples, not an exhaustive search.

These results are finite exact checks, not a proof of the universal fixed-dimension lower bound or Lean kernel verification. The proposed general proof uses test-vector-dependent merging, Serre's strict-halfspace lemma and the published Lai–Lim count bound. No Lean proof is included.
