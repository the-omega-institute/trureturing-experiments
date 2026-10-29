# Javaheri cyclic arrangement dimension

Exact arithmetic checks for Mohammad Javaheri, *A new arrangement inequality
(revised and extended version)*, arXiv:2610.33653v1, Conjecture 1.3.

Run with Python 3 and SymPy (checked with SymPy 1.14.0):

```sh
python3 docs/reports/javaheri-cyclic-arrangement-dimension/check.py
```

The stored output is `check.out`; the final line is `ALL_CHECKS_OK`.

The source fidelity check reproduces Corollary 4.2(ii), label `cor:examples`,
at b=3: N(0,3)=1. The source inequality is equality in dimension 1, while
(2,1,...,1) has left side minus right side 2(n-1) for every n>=2.

For the proposed refutation, the program verifies the two-variable
factorization and its nonnegative-square decomposition at (a,b)=(6,1), then
verifies the all-dimension identity

    D_n(6,5,...,5) = 91(1025n-3071) = 93275(n-3)+364.

Thus the displayed witness is positive for every n>=3. Direct evaluation in
dimensions 3 through 100 checks cyclic indexing independently of the symbolic
sum formula. At n=3 the cyclic sum is 405155, the left side is 6482480,
the right side is 6482116, and the difference is 364. The program also checks
the family (m+1,m,m) for integers m>=5, including a polynomial identity with
strictly positive coefficients after substituting m=r+5.

The exact threshold values 48 sin(pi/n)^2 at n=2,3,4 are 48,36,24,
respectively, against the required value 36. The global exclusion of n>=4
uses monotonicity of sine on [0,pi/2]: sin(pi/n)^2<=1/2. This analytic step
is stated explicitly rather than inferred from a finite numerical sample.

The normalized family (1+t,1,1) has zero quadratic coefficient and positive
cubic coefficient in D_3. It explains the obstruction at the local
second-variation boundary. The conjecture asks for global validity.

All computational evidence uses integer, rational, or exact symbolic
arithmetic. There is no floating point evidence and no Lean proof here.
The proposed mathematical conclusion is N(6,1)=2, whereas N*(6,1)=3.

Source: https://arxiv.org/abs/2610.33653v1
