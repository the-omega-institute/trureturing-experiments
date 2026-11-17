# Fixed-parent FIB exceptional layers

The coefficient and norm identities come from `fib_shared_response_177/verify.py`,
with the original bundle `fib_shared_response_176/global_preflight` at commit
`7c45a21dbea5cdc0974f878a96dc64281cc5db07`. The supplier's source identities,
literal clock convention, complete cubic gain recovery and image resultant are
checked before the exceptional-layer calculations.

Run with Python 3.12 or later, SymPy 1.14.0 and python-flint 0.9.0:

```sh
python3 fib_shared_response_178/verify.py
```

The certificate gives the exact valuations of all six remaining factors in the
source and image discriminants and chart resultants. Factors of degree 4, 24 and
32, with 48 real roots in total, do not divide the image discriminant. This fact
requires a separate original-source coverage argument before it implies response
separation at a chart exception.

For each of the three quadratic factors the actual reduced source has degree 14.
The first uses the actual leading-degree drop; the second removes only the
factor `(u+1)^2`; the third removes only `u^2`. The latter two factors cannot
represent a legal positive source. Exact coefficient reduction in the quadratic
field and degree-preserving finite-field specialization give units for A, D and
Phi-sharp, a squarefree reduced source, and a cyclic MM-response algebra. The
cyclic determinants are 98 modulo 101 at b=59, 54 modulo 101 at b=10, and 58
modulo 103 at b=44. The full matrices and Bezout identities are in the certificate.

These are algebraic certificates for the fixed parent theta=3/7, xi=-35/29.
They do not decide the original adaptive root problem and are not Lean proofs.
