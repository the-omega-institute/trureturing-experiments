# Baluyot–Conrey–Irfan–Snaith Conjecture 1: exact preregistration checks

Source: arXiv:2611.05000v1, submitted 3 November 2026, Conjecture 1,
equation (1.16), and the polynomial range in Theorem 2.2.

Run `python3 docs/reports/baluyot2611c1/check.py` with Python 3.12 or later
and SymPy. The stored standard output is `check.out`; success means exit 0
and the final line `ALL_CHECKS_OK`.

The fidelity check evaluates the source's defining residues for the two
polynomials P_{0,1} and P_{1,1}, and recovers its stated piecewise expression
(1.15) through Theorem 2.2. The independent check evaluates the k=1 moment
contour exactly, checks the symbolic main-term gap for averaging heights
N^q and N^(q+1), verifies admissibility in 594 finite batches, certifies 33
strictly positive gap coefficients by rational Taylor enclosures, and checks
the exact finite-average error identity and its sharp factor 1/2 for sizes
2 through 30. All calculations use integers, rational numbers, or symbolic
expressions. No floating point arithmetic is used as evidence.

The obstruction is general: if B is the average of N translated short
averages A_r, whereas their proposed main terms are Q and P respectively,
then Q-P = average(A_r-P) - (B-Q). Thus at least one error is at least
|Q-P|/2. Here Q-P is a positive constant times (log N)^2, while the
conjecture allows only O((log N)^(1+xi)), with xi<1.

The script does not evaluate zeta, locate a zero-free line, integrate zeta,
or prove a Lean theorem. Its rational choices of the real parameter check
algebra and range constraints only. The proposed unconditional mathematical
family chooses the parameter outside a countable exceptional set of zero
abscissae; that analytic existence argument, contour evaluation in Lean,
and the asymptotic contradiction await the post-filing probe. No claim or
negation has been proved in Lean by this preregistration seat.
