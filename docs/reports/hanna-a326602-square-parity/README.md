# Hanna A326602 square parity

Run `python3 check.py` with Python 3. It uses only standard-library exact integer
arithmetic and requires no network access. Expected output is in `check.out`.

Source: [Paul D. Hanna, OEIS A326602](https://oeis.org/A326602), revision 19,
FORMULA and the 68 displayed terms. The source's coefficient convention is
`A(x) = sum a(m) x^(2m)`; its index-zero summand is `1/(1+x)`.

The program expands FORMULA (1) through degree 800 by exact polynomial inversion,
then independently enumerates the finite triples with exponent
`l+k+2(l+j)(j+k)` and multinomial weight. It compares all coefficients, checks the
swap of `l` and `k`, checks the fixed-point central-binomial factor, and verifies
that all 400 odd-degree coefficients vanish and that `a(0)` through `a(400)` are
odd exactly at square indices. All 68 source terms are reproduced.

These finite checks support the proposed unbounded argument. They are neither a
proof of the conjecture nor a Lean probe or a formal settlement.
