# Rational digits and simultaneous convergence

The check targets Conjecture 16 in Giuliano Romeo, *Real convergence and periodicity of p-adic continued fractions*, arXiv:2410.09215v3, and The Ramanujan Journal 68, article 112, DOI 10.1007/s11139-025-01264-7.

Run `python3 check.py` with Python 3. It uses only integers and `fractions.Fraction` from the standard library. The stored output is `output.txt`.

For source fidelity, the program checks the source's continuant determinant identity on every word of lengths 2 through 5 over the digits −2, −1, 0, 1, 2 and 1/2. It also checks the source's golden mean quadratic fixed-point identity.

For the proposed refutation, exact arithmetic in Q(phi), where phi² = phi + 1, constructs two nine-digit prefixes at p = 11 and one at p = 19. The p-adic embeddings use the simple roots phi = 4 modulo 11 and phi = 5 modulo 19. Real comparisons use certified rational isolating intervals; p-adic valuations use Hensel lifts with a precision criterion certifying each returned valuation. The first p = 11 prefix starts with 14/9.

At every checked step n, the program verifies the complete-quotient selection conditions, the exact convergent error identity, positivity, the factorial denominator bound (n ≥ 1), the p-adic denominator bound, and error valuation at least (n + 1)². The selected digits satisfy a_n > n. Distinct selection rules produce distinct prefixes.

These are finite exact checks of the proposed construction and its identities. Infinite existence, convergence and failure of eventual periodicity require the general argument stated in the preregistration. The construction permits arbitrary rational digits; it imposes no Browkin digit legality or restriction to Z[1/p].
