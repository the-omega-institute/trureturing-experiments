# Alternating-chord exterior path

Source: Pavel Chebotarev, *How the eigenvalue region of standardized Laplacians leaves the 2005 polygon*, [arXiv:2611.02249v1](https://arxiv.org/abs/2611.02249v1), Section 5, the numerical assertions after Proposition 10.

The target has even integers n,m with m >= 2 and 2(m+1) < n. It concerns the continuous root of the stated two-dimensional restriction starting at lambda_m, its endpoint q_m, and strict exteriority to S_n for 0 < w < 1.

## Exact check

Requires Python 3 and SymPy 1.14.0. From this directory run:

```sh
python3 check.py > check.out
```

The program uses rational coefficients, exact symbolic algebra and reduction modulo cyclotomic polynomials. It uses no floating point. It checks the quadratic characteristic polynomial, its discriminant, the crossing-exclusion identity, the second-order departure coefficient and the algebraic decomposition used for Re R >= 3/4. For all 81 admissible pairs with even n <= 40 it reproduces both columns of Proposition 10's Fourier restriction directly from the original matrix row formula, the paired endpoint eigenvector, both endpoint roots, the trace and the normalized discriminant. This includes 2,352 original-matrix rows. The final line is `ALL_CHECKS_OK`; exit code 0.

The infinite-family trigonometric inequalities and continuation argument are mathematical obligations in the preregistration, rather than conclusions certified by this program.

## High-precision diagnostic

Requires mpmath. Run:

```sh
python3 check_high_precision.py > check_high_precision.out
```

At 120 decimal digits, the diagnostic evaluates the explicitly selected square-root branch at w=j/1000 for j=0,...,1000, for every admissible pair with even n <= 40. It checks the characteristic residual, initial and terminal roots, Re R >= 3/4, strict exterior supporting-line sign at all 80,919 interior samples, and the inward sign at every polygon vertex. The final line is `HIGH_PRECISION_CHECKS_OK`; exit code 0. These finite numerical readings support the proposed argument and do not certify the continuum or the unbounded family.

The source's numerical lower bound 0.14 on root separation is not a target of either program and is not asserted uniformly in n.
