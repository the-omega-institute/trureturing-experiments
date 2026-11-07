# Abed–Nikolaeva–Novikov, Conjecture 11.1, item 9

Source: Sami Abdullah Abed, Irina Nikolaeva and Andrej Novikov, *How characterisations of traces and central elements fail: sharp violation constants via Michelson contrast and noncentrality*, [arXiv:2610.34388v1](https://arxiv.org/abs/2610.34388v1), Proposition 5.17 and the item-9 clause of Conjecture 11.1.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The stored output is `output.txt`. Successful execution ends with `ALL_CHECKS_OK` and exits 0. All evidence uses exact integer, rational, complex algebraic and symbolic arithmetic; no floating point is used.

For positive definite complex matrices rho and self-adjoint idempotents P,Q, define F = Tr(rho(PQ+QP−2QPQ)) and H = Tr(rho(P+Q)). The proposed envelope at contrast delta is s/(1+s), with s² = (1−sqrt(1−delta²))/2. No trace-one normalization is imposed; H=0 is excluded.

The check reproduces Proposition 5.17's coordinate identity, its slope 1/2 at zero and its limiting value sqrt(2)−1. It verifies the symbolic equality family R = [[t,2s sqrt(t)],[2s sqrt(t),2−t]], P = [[1,sqrt(t)],[sqrt(t),t]]/(1+t), Q = diag(1,0), with t=1−2s²>0 and 0≤s<1/sqrt(2). The characteristic polynomial gives eigenvalues 1±delta, where delta²=4s²(1−s²); the quotient is exactly s/(1+s). At delta=24/25 it checks F=21/40, H=7/5 and F/H=3/8, including zero padding through dimensions 2–7 and the delta=0 endpoint.

The finite regression sets contain 1,584 admissible qubit cases and 999 admissible rank-mask cases in dimensions 2–4, with all ranks represented. Their finite sets of projection frames and contrasts are explicit in the script; 25 pairs with P=Q=0 are excluded and checked separately. A further 27 cases verify trace additivity and contrast contraction for compressions of the same density, including densities with nonzero off-diagonal blocks. Every inequality is compared by rational arithmetic with sign guards before squaring.

These checks support a preregistration with expected settlement **Proved**. They do not prove the all-dimensional upper bound or constitute a Lean proof. Other clauses of Conjecture 11.1 and the singular boundary delta=1 are outside scope.
