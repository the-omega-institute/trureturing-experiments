# Hanna A352409 parity checks

Run `python3 check.py` with Python 3. The program uses only integers and
`fractions.Fraction`, requires no packages or network, and exits with
`ALL_CHECKS_OK` when every assertion passes. Its captured output is `check.out`.

The source is Paul D. Hanna's [OEIS A352409](https://oeis.org/A352409),
submitted March 15, 2022, revision 17 (April 23, 2026), NAME and COMMENTS (4).
The embedded 16-term prefix and coefficient table are attributed to Hanna;
OEIS content is subject to the [OEIS End-User License](https://oeis.org/LICENSE).
The companion is [A352408](https://oeis.org/A352408), revision 12.

The program reproduces the displayed prefix and seven coefficient-table rows.
It constructs 40 coefficients by the defining triangular equations, cross-checks
power coefficients using two independent algorithms, and checks oddness and
the defining zero diagonal. For indices 1 through 100 it compares a direct
power calculation with the finite reference-residual expansion, checks the
leading binomial identity, and checks evenness term by term in the rationals
with odd denominator. It also tests normalized congruence transfer on 188
pairs with fixed seed 352409.

Modulo-three checks reproduce six companion terms, test its residuals through
100, and compare the 40 generated coefficients with the Catalan residue
pattern. These are finite observations: the uniform modulo-three residual
identity remains unproved. No check here proves a uniform assertion or
constitutes Lean verification. The preregistration target is A352409 COMMENTS
(4) only; the planned settlement uses a general argument, with these checks
serving as independent arithmetic and source-fidelity checks.
