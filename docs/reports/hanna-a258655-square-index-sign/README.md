# Hanna A258655 square-index signs

Run `python3 check.py` in this directory. Python 3's standard library suffices;
all evidence uses integers and exact rational arithmetic. The output is in
`output.txt` and ends with `ALL_CHECKS_OK` on success.

The source is Paul D. Hanna's [OEIS A258655](https://oeis.org/A258655/internal),
revision 20, with its parent [A256357](https://oeis.org/A256357/internal),
revision 26. The included `b258655.txt` is the source's 1000-term
[table](https://oeis.org/A258655/b258655.txt), SHA-256
`f9638076cc38cc9467a77ce6bbd5b7e19db330dd4e6e7a6afc67566560fa0c89`.
Source attribution and reuse follow the [OEIS End-User License](https://oeis.org/LICENSE).

The checks reproduce all table values, compute the NAME's logarithmic
coefficients independently through degree 10000 using the triangular identity
`q F' = F sum(B_N q^N)`, and cross-check that recurrence with the defining
rational logarithm through degree 64. They compare the published divisor
formula, odd-character factorization, two-adic decomposition and strict signs.
They cover all square indices n = 1..100 from the NAME and all n = 1..1000
from the divisor formula, including the odd part 1 and two-adic exponents 0..8.
An additional exact convolution check verifies through degree 512 that
`C(q)^2 = F(q) psi(-q)`, where `C(q) = sum over k in Z of (-1)^k q^(4k^2-k)`,
and that `psi(q) phi(-q^2) = psi(-q) phi(q)`. These are the proposed bridge's
theta identities, independent of the divisor formula.

These are finite checks. They neither prove the uniform conjecture nor supply
a Lean proof of the theta-product identity connecting the NAME to the divisor
formula. The planned settlement is a general character and divisor argument;
these computations are supporting evidence only.
