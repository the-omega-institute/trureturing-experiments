# Modified Bessel ratio unimodality: exact checks

Target: the first sentence of Conjecture 1 in D. Karp, A. Vishnyakova and
Y. Zhang, *Unimodality preservation by ratios of functional series and
integral transforms*, arXiv:2408.01755v2, Section 4, printed page 12,
DOI 10.1007/s00025-025-02431-4. The target allows real orders
nu1 >= nu2 > -1 and scales 0 < a1 <= a2. The separate log-concavity
clause is outside this report.

Run with Python 3 and SymPy 1.14.0:

```sh
python3 check.py
```

The output is in [check.out](check.out). All evidence uses integer and
rational arithmetic or exact symbolic algebra; no floating-point evaluation
contributes to a check. Failures raise an exception even under Python's
optimization flag.

The source-fidelity checks reproduce Section 4's equal-order monotonicity:
normalized coefficient ratios are (a1/a2)^(2k), and exact intervals certify
nonpositive quotient derivatives at the specified positive points.

The remaining checks verify the modified-Bessel coefficient recurrence,
the hypergeometric derivative and contiguous coefficients, the weighted
Wronskian identity, the quotient derivative identity and the leading endpoint
coefficients for equal and unequal orders. The grid includes negative orders,
equal orders, equal scales and noninteger order gaps. It checks all 36 ordered
pairs from the eight listed orders, four scale pairs and eight positive points.
The 128 literal constant cases are distinguished from interval sign decisions.
Every nonconstant sign decision encloses the infinite series, using a geometric
tail after 128 terms. The derivative signs have no negative-to-positive
transition on this grid.

For b > 0 and t > 0, the interval routine encloses
F_b(t) = sum t^k/(k! (b)_k). Its omitted first term is t_N and its subsequent
term ratios are at most rho = t/((N+1)(b+N)) < 1, giving the upper tail
bound t_N/(1-rho). The derivative identity
x I_nu(a*x)'/I_nu(a*x) = nu + 2t F_(nu+2)(t)/((nu+1)F_(nu+1)(t)),
with t = a^2*x^2/4, turns these intervals into derivative-sign certificates.

These are finite diagnostics and exact algebraic identities. They do not prove
the universal conjecture or provide a kernel-checked settlement. The proposed
general argument uses the sign changes of the weighted Wronskian and an
endpoint argument. The source's displayed ratio expansion contains the
prefactor 2^(nu2-nu2); the standard defining series instead gives
2^(nu2-nu1). The checks use the standard series, as does the repository's
existing Bessel definition; they do not adopt that typographical prefactor.
