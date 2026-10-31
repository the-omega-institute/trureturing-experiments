# Yang oscillating series: exact finite checks

Source: Kai Yang, *Analysis of Convergence for Some Positive Oscillating Series and Integrals*, Memorial University MSc thesis, June 2011, Lemma 4.1 and Section 6.1 Question 5. Institutional record: <https://hdl.handle.net/20.500.14783/1902>.

Run `python3 check.py`. Python's standard library is sufficient. `output.txt` is its output. All comparisons use integers or `fractions.Fraction`; no floating point calculation is evidence.

The fidelity check reproduces Lemma 4.1 for theta equal to the square root of 2 and every integer N from 1 to 512, deciding the inequalities by exact squaring. The other checks exercise the proposed proof's absolute-value kernel bound, perturbed rational orbit counts (including closed endpoints and wraparound), and six dyadic block mass bounds for the quadratic phase alpha equal to the square root of 2 minus 1. The latter uses an integer-certified rational enclosure of the phase and rounds rational upper bounds upwards.

The block mass calculation uses the rational majorant `1/(1+8*kappa*M*d^2)` that follows from the proposed Gaussian estimate, with a equal to 1/4, p equal to 1 and kappa equal to 1/5. It checks this majorant directly; the transcendental Gaussian estimate remains a mathematical proof obligation. None of these finite checks proves the infinite convergence assertion, the Diophantine membership of the example, or absence of prior art. Positive settlement requires a general argument, with these computations serving only as checks of helper statements.
