# Kimberling A389142 urn-likelihood checks

Run `python3 check.py` with Python 3 and SymPy. The program uses integers,
`fractions.Fraction`, and exact symbolic polynomials throughout. Its saved
standard output is `check.out`; a successful run ends with `ALL_CHECKS_OK`.

The source is Clark Kimberling, [OEIS A389142](https://oeis.org/A389142),
revision 11 (2025-10-05), and Marc A. Brodie,
[OEIS A111097](https://oeis.org/A111097), revision 35 (2026-04-02).
Source numerical data are attributed to the OEIS Foundation and contributors;
see the [OEIS End-User License Agreement](https://oeis.org/LICENSE).

The program reproduces all 54 displayed A389142 terms and the A111097 example
at n=3. It checks the finite-draw counting interpretation, unique maxima for
n=2..60 over k=1..(2M+10), the two boundary signs through n=300, exact
binomial truncation bounds, six polynomial identities, and the likelihood
ratio and geometric identity. The maximum scan has a finite cutoff; it
makes no claim about larger k or an all-parameter Lean proof. The symbolic
identities alone do not prove the statistical maximum.

The proposed settlement requires a general argument for every n>=2 and every
positive integer k. The compiled preregistration conventions contain only
definitions and the proposition, with no proof of it or its negation.
