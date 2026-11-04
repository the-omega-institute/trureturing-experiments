# Dimitrov: reciprocal-product bound for φ, ψ and σ

Source: S. I. Dimitrov, *Lower bounds on expressions depending on the functions φ(n), ψ(n) and σ(n), III*, [arXiv:2606.12484v1](https://arxiv.org/abs/2606.12484v1), Conjecture 1 and Theorem 1.

Run with Python 3 and SymPy (checked with SymPy 1.14.0):

```sh
python3 docs/reports/dimitrov-phi-psi-sigma-reciprocal-bound/check.py
```

The program uses integers, exact rational fractions and symbolic rational algebra. It reproduces Theorem 1 for every integer from 2 through 300, verifies the conjectured inequality and the proposed nonnegative-gap decomposition for all 3,588 pairs with 2 ≤ n ≤ 300 and 1 ≤ k ≤ 12, and checks 676 independent rational triples including the R=1 boundary. It checks the source's arithmetic-function products against counting definitions and the existing D5 representation of Dedekind psi. The symbolic residual of the proposed gap identity is exactly zero.

The output is in `check.out` and ends with `ALL_CHECKS_OK`. These checks support a proposed symbolic proof; they do not prove the universal assertion in Lean or establish literature novelty. The finite search is corroboration, and is not the proposed positive settlement's main content.
