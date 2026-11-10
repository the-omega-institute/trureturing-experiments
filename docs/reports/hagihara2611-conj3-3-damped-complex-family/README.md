# Hagihara Conjecture 3.3: a damped complex family

Source: Fugen Hagihara, *Positivity Proofs for Third-Order Holonomic Sequences via Differences Scaled by Rational Functions*, arXiv:2611.10920v1, Conjecture 3.3 (PDF; Conjecture 7 in HTML).

`check.py` uses exact integer, rational and SymPy algebra. Run with Python 3.12 or newer and SymPy (`python3 check.py`). Its standard output is stored in `output.txt`; success ends in `ALL_CHECKS_OK`. The check was run with Python 3.14.8 and SymPy 1.14.0.

The fidelity check reproduces the Franel-sum recurrence stated in source Example 3.2 at indices 1 through 30. The settlement checks verify the common recurrence of the family with `4/3 < b^2 <= 5`, its initial-value determinant, a quantitative positive cone at `b=2`, and the failure of the conjectured necessary coefficient bound. Three product solutions are checked independently against the recurrence at small indices.

For `b^2=2`, the original recurrence is iterated directly through index 10,000 from two initial triples: `(1,4/5,1/2) + (1,-1,1)/2048` and `(1,4/5,1/2) - (1,-1,1)/2048`. Every computed term is strictly positive. These iterations use `fractions.Fraction` and no closed-form products or floating-point arithmetic. Output includes trace digests and bit sizes instead of unwieldy full terms.

These finite checks support a prospective refutation; they do not prove the infinite statement or certify a Lean theorem. The infinite argument uses a factorwise telescoping bound and an invertible initial-value map. The source's first-order-product construction is an acknowledged ingredient. The limiting characteristic polynomial has a triple root; the complex product modes decay relative to a positive solution.
