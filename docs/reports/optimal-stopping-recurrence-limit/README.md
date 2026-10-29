# Pointwise recurrence-to-ODE convergence

Exact checks for Conjecture 1, Section 5.1, of L. Bayón, P. Fortuny Ayuso,
J. M. Grau, A. M. Oller-Marcén and M. M. Ruiz, *A new method for computing
asymptotic results in optimal stopping problems*, arXiv:2205.08495v1.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The captured result is
in `output.txt`. All numerical evidence uses integers and `fractions.Fraction`;
all symbolic checks use exact SymPy expressions.

Fidelity is checked independently by enumerating every permutation of up to
eight applicants, calculating each threshold's success probability, and testing
the source's classical-secretary recurrence. The source's limiting ODE solution
`-x log(x)` is checked symbolically.

The candidate counterfamily uses parameters `0 < alpha, theta, q < 1` and
`0 < c < (1-theta)q`. A spike lies at
`j_n = floor(alpha*n) + isqrt(n)` when `1 <= j_n < n`.
For `k < n`, put `B=(n-k)/(n-k+1)`,
`H=B*(1-q*indicator(k=j_n))` and
`G=theta*(1-H)+c*indicator(k=j_n)`, with terminal value `theta`.
The script solves each row backwards and independently checks its closed form,
probability bounds, residual concentration, and exact rational interior values.
It also checks the ODE invariant and the predicted jump symbolically.

The intended analytic limit is `theta+c*(1-alpha)/(1-x)` for `0<x<=alpha`
and `theta` for `alpha<x<1`. The coefficient limits are `theta/(1-x)` and
`1/(1-x)`, whose ODE has only the constant finite-endpoint solution `theta`.
Finite row checks and symbolic substitutions do not prove the universal limit,
ODE classification, or the conjecture's negation. No Lean proof is supplied.
