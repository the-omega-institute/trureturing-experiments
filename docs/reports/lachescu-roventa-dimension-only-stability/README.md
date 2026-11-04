# Dimension-only stability of omega-m-star-convexity

Exact-arithmetic supporting checks for Lăchescu–Rovenţa, arXiv:2207.08058v1,
Section 5, Problem 1. The source asks for an exact approximation with a constant
depending only on dimension. The proposed negative answer uses C=[0,1],
m=1-1/(2n), omega(t)=sqrt(t), delta=m/4, and Phi(x)=x. The exact comparison
class contains Psi_0(x)=x-m/(1-m). The accompanying general argument forces
k_1 >= sqrt(n/2) for every n>=1.

Run with Python 3 and SymPy:

```sh
python3 docs/reports/lachescu-roventa-dimension-only-stability/check.py
```

`output.txt` is the program's output. Arithmetic uses integers, Fraction, and
symbolic algebra; radical inequalities are compared through positive integer
powers. No floating-point calculation is evidence.

Fidelity checks reproduce the source's barycentric proof coefficients for
convex functions with f(0)<=0, the square-function gap certificate, and the
nonconvex polynomial example's midpoint gap 5/256. Further checks cover the
approximate/exact gap identities, special-pair cancellation, symbolic telescopes
through length 16, and the lower-bound arithmetic for n<=64 and four rational
exponents q in (0,1). The finite grids are regression evidence, not proofs over
all real inputs. The general refutation requires the ordinary argument and a
subsequent kernel-checked Lean proof; no Lean settlement proof is supplied here.

Source: https://arxiv.org/abs/2207.08058v1 (submitted 17 July 2022),
DOI https://doi.org/10.1007/s00010-023-00942-5.
