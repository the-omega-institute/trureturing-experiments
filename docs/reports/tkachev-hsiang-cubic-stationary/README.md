# Stationary harmonic cubics outside the Hsiang class

Source: Vladimir G. Tkachev, *Hsiang cubics as spherical Chebyshev polynomials:
equioscillation and idempotent quadrature*, arXiv:2610.13532v1,
Problem 11.3(b), <https://arxiv.org/abs/2610.13532v1>.

Run `python3 check.py`. Python 3 and SymPy are required. The recorded output is
`output.txt`; a successful run ends with `ALL_CHECKS_OK`. Every algebraic value
and comparison uses exact integers, rational numbers or algebraic expressions.

The source-fidelity check reproduces the cubic of Pol(1,4), its harmonic Hsiang
identity, squared spherical norm 2/315 and Chebyshev ratio 35/6 (source labels
pol14 and tab:ratio), and the monomial normalization in Lemma 6.1.

For m = 2, 3, 4 the program checks the orthogonal sum of two copies of
q(t,y,z) = t(sum y_j^2 - sum z_j^2) in dimension 4m+2. It verifies harmonicity,
the full infinity-Laplacian polynomial, the 8m positive extremal nodes, their
complete cubic moment polynomial, and the cubature/Fischer/spherical identities
against every ambient cubic monomial, including mixed-block directions. It
also checks a rational parameterization of infinitely many maximum points,
two evaluations incompatible with every Hsiang parameter, and a trivial kernel
for the map taking a constant vector to its directional derivative of p.

The nonnegative factorization underlying the exact supremum is checked
symbolically. Infinitude, the one-sided directional-derivative limit defining
stationarity, and the family for every m >= 2 require the general mathematical
argument in the preregistration; finite checks do not prove those assertions.
There is no Lean proof in this report.
