# Two separable qutrit probes

Exact algebraic checks for the unnumbered closing question in Section VII of P. Ghosal, P. Halder, A. Patra and A. Sen(De), *Quantum hypothesis testing of non-mixed-unitarity: A multifaceted hierarchy of quantum channel discrimination*, [arXiv:2609.40355v2](https://arxiv.org/abs/2609.40355v2).

Run with Python 3 and SymPy 1.14.0:

```sh
python3 check.py
```

The program uses exact integers, rationals, radicals and symbolic polynomials. It reproduces the source's Theorem 6 output and compression for the antisymmetric two-qutrit probe. It verifies an explicit eleven-product-projector decomposition of `(I9 + 2 rho)/11`, the whole symbolic dilution family, its partial-transpose spectrum, and the finite algebraic certificates used by the compression obstruction. The dilution identity is checked with 72 independent symbolic coefficients describing an arbitrary unital complex-linear map.

The stored output is in [output.txt](output.txt), ending in `ALL_CHECKS_OK`. These checks support an affirmative proof outline. They do not enumerate all mixed-unitary channels, compute a positive numerical separation constant, prove compactness or Pinsker's inequality, certify the operational Stein exponent, or provide a Lean proof.
