# Exact checks for Sra's Koteljanskii questions

Source: Suvrit Sra, *A Koteljanskii inequality for permanents*,
[arXiv:2609.39979v1](https://arxiv.org/abs/2609.39979v1), Section 4(b,d).
The matrix in the fidelity test is the published Theorem 1.2 matrix.

Run with Python 3 (standard library only):

```sh
python3 docs/reports/sra-koteljanskii/check.py
```

`check.out` is the exact output. Successful execution ends in `ALL_CHECKS_OK`.
The program uses integer permutation sums and rational matrix arithmetic.

- Fidelity: verifies the source's inverse, four permanents, and product ratio
  `832/833` in Theorem 1.2.
- Item (d): constructs the whole reciprocal-permanent polynomial, differentiates
  in the three common indices, and verifies a positive definite restriction of
  its Hessian to a two-dimensional plane. Such a restriction excludes the
  necessary Lorentzian signature condition. This is a new application of the
  published matrix, rather than a new discovery of that matrix.
- Item (b): verifies a symmetric positive integer matrix and its inverse M-matrix
  certificate. Exact cycle enumeration gives the supermodularity difference
  `alpha^2*(alpha+1)*(2401*alpha^2+1617*alpha-16)`.
  Its bracket is at most `-314199/40000` throughout `0<alpha<=1/200`, because
  its nonconstant coefficients are positive. Thus the parameter interval is
  certified by a monotonic polynomial bound, rather than by sampling.
  The program also checks all subset pairs at the known orders `1/2` and `1`.

These are exact computational certificates and mathematical proof outlines;
they contain no Lean proof of either settlement. Item (b) disproves the universal
extension; it does not classify the remaining positive orders. The Lorentzian
necessary condition is from Brändén–Huh, *Lorentzian polynomials*, Definition 2.1
and Example 2.4, [arXiv:1902.03719v8](https://arxiv.org/abs/1902.03719v8),
DOI [10.4007/annals.2020.192.3.4](https://doi.org/10.4007/annals.2020.192.3.4).
