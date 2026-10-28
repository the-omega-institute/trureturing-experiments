# Rising-factorial Hurwitz stability

`check.py` checks the proposed refutation of D. Karp, *Positivity of
Toeplitz determinants formed by rising factorial series and properties of
related polynomials*, arXiv:1203.1482v1, Conjecture 2, equation (3).

Run with Python 3 and SymPy:

```sh
python3 docs/reports/rising-factorial-hurwitz/check.py
```

`output.txt` is the complete successful output. All evidence uses integers,
`fractions.Fraction`, and exact SymPy algebra. There is no numerical root
finder or floating-point evidence.

The fidelity checks reproduce Section 2, Lemma 3 (`f_k=1` gives the zero
polynomial) for symbolic shifts and degrees 0 through 8, and independently
compare the determinant-series coefficient with equation (3) through index 5.

For `n=2m`, shifts `alpha=beta=A`, and `f_k=A^(-(k-m)^2)`, the program
constructs the integer polynomial `A^(2m^2) Q_(2m)^(A,A)`. It checks strict
log-concavity and positivity of every coefficient. Exact Routh sign counts
give two right-half-plane zeros for `(m,A)=(10,1000),(10,10000),(11,1000),
(12,1000)`. Routh rows are divided only by positive factors; zero pivots or
all-zero rows cause failure. Four polynomials with explicit elementary roots
check the sign-count convention. At `(m,A)=(10,100)`, the right-half-plane
count is zero, recording a finite parameter control for the necessary
large-parameter qualification.

The limiting polynomial `(z+1)^(2m)-(z(z+2))^m` has degree `2m-2`, leading
coefficient `m`, and positive coefficients for the checked range `2<=m<=16`.
Its unstable cyclotomic root at `m=10` is checked algebraically, including
the exact inequality `sin(pi/10)=(sqrt(5)-1)/4<1/3`. The finite instances and
identities support a proposed asymptotic family; they do not constitute
kernel verification or a check of every real parameter. Uniform convergence
and persistence of the unstable zeros remain obligations for formalization.
