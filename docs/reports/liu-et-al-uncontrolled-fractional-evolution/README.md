# Uncontrolled fractional evolution: exact checks

Source: Liu, Wesołowski, Thompson, Gu and Pira, *Optimal query complexity for
fractional quantum evolution*, arXiv:2610.01940v1, Section V, “Uncontrolled access”.

Run with Python 3 and SymPy 1.14.0:

```sh
python3 check.py
```

The program reproduces the Laurent support bounds in the proof of the source's
scalar-polynomial reduction on all 15 query-sign words of length at most three.
The hard family uses the exact admissible gap pi/4 and nontrivial inter-query gates.

For the proposed negative answer, the program checks the nonscalar oracle pair
U = diag(i, -i), V = -U, its prescribed Hermitian logarithms, and symbolic target
outputs for arbitrary real t. Their difference has squared singular values
sin(pi t)^2, so its trace norm is 2|sin(pi t)|. On 0 < t < 1 this becomes
2 sin(pi t). The exact t = 1/2 outputs are orthogonal pure states; at epsilon = 1/4
their separation 2 exceeds twice the permitted error, 1/2. It also checks 62
circuits with both query signs, zero through four queries, and workspace
dimensions one and two. At t = 0 and t = 1 the target-channel separation vanishes.

The spectral and gap comparisons use 3 < pi < 4, certified by the exact
eight-term Gregory lower partial sum for pi/4 and its first-term upper bound.
The checker uses only rational comparisons for these inequalities.

All calculations use exact integers, rationals, algebraic constants and symbolic
expressions. The stored output is `check.out`, ending in `ALL_CHECKS_OK`.
Finite enumeration does not establish the arbitrary-query statement. That step
uses induction on circuit composition; the final error inference uses the
unnormalized diamond norm and the trace-norm triangle inequality. These remain
proof obligations for the subsequent Lean probe. No Lean proof is included.
