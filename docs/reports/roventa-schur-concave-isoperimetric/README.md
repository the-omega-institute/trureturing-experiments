# Rovenţa’s Schur-concave isoperimetric question

Source: Ionel Rovenţa, *A note on Schur-concave functions*, Journal of
Inequalities and Applications **2012**, 159, Section 2, equation (2.3) and
its final open question. DOI: https://doi.org/10.1186/1029-242X-2012-159.
The article is licensed CC BY 2.0. The check is an independent implementation;
it contains no copied article source.

Run `python3 check.py` using Python 3's standard library. The stored result is
`output.txt`; successful execution ends with `ALL_CHECKS_OK`.

For fidelity, the program reproduces Theorem 5 for the positive concave function
f(t)=t/(1+t) on 1,089 inputs, and the exact additive identity in (2.6). Ratios
are positive, so the rational comparison A^(n−1) ≤ B^n is equivalent to (2.3).

For the expected refutation it uses F=M_(−p), with
M_(−p)(y)=((sum y_i^(−p))/n)^(−1/p). It checks x=(1,…,1,n+1) for
3 ≤ n ≤ 16 and p=n² or 2n², six inputs with all coordinates distinct,
and three other inputs having a repeated minimum. Every root is enclosed by
rational bisection with certified p-th-power comparisons. The program uses no
floating point. It also checks finite T-transform inequalities, the strict
Schur–Ostrowski signs after factoring a positive common coefficient, and an
exact weighted Cauchy–Schwarz identity controlling the Hessian.

These are finite checks and algebraic consistency checks. They do not prove
Schur-concavity on the full positive orthant, the universal family in n and
real p, or novelty in the literature. The intended universal argument is the
minimum envelope a ≤ M_(−p)(y) ≤ n^(1/p)a: when at least two coordinates
attain a, every replacement retains a, and exponent separation forces failure
if (mean(x)/a)^p > n^n. No Lean proof is included.
