# Three-state Bargmann density: exact supporting checks

Source: Lin Zhang, *Exact probability laws of random inner products: Dirichlet
splines, overlaps, uncertainty, and Bargmann invariants*, arXiv:2610.25196v1,
Section 6 and the second open-problem bullet in Section 7.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The program uses Python
integers, `fractions.Fraction`, and exact SymPy algebra, integration, and limits.
It does not use floating point or numerical quadrature. The result is stored in
[output.txt](output.txt); success ends with `ALL_CHECKS_OK` and exit code 0.

The fidelity check reproduces the source's qubit first and second marginal
moments. The proposed disk-density reduction is checked against all 175 mixed
moments for dimensions 3 through 9 and exponents 0 through 4. Additional exact
checks cover normalization, logarithmic triangle bounds for exponents 0 through
8, the boundary Hessian and Jacobian, quadratic integrals for dimensions 3 through
11, and 79 rational boundary parameters, including both real endpoints.

These checks support the proposed general proof. They do not certify Haar
conditional laws, changes of variables for singular integrals, continuity of the
density representative, coercive sublevel localization, or the full asymptotic
statements. The intended settlement uses an argument in every integer dimension
at least 3; the finite checks are supporting evidence.
