# Niu Question 13: exact preregistration checks

Source: Tong Niu, *The first descent in a standard Young tableau of shape (n,n,n)*,
[arXiv:2610.03038v1](https://arxiv.org/abs/2610.03038v1), Question 13.

Run `python3 check.py` with Python 3 and SymPy. The stored reading is
`check-output.txt`. The program uses integers, fractions and exact symbolic algebra;
it has no network or source-file dependency.

The program enumerates all standard tableaux for n = 1,...,5, checks the source's
initial values, and compares the recurrence and finite-sum formula through n = 140.
It prints the exact recurrence values through n = 40. Full column rank over Q
excludes complex-coefficient first-order fits of degrees 0,...,12 on n = 2,...,39
and n = 20,...,100. These finite windows do not exclude an eventual relation
starting after the windows, or one of larger degree.

The remaining checks verify the proposed alternating hypergeometric solution,
the order reduction, the Gosper factorization, all nine rational root gaps,
the leading coefficient obstruction, the source telescoping certificate, and
the corner recurrence used for the general three-row cardinality formula.
These are checks for a proposed general argument. They are neither a Lean proof
nor a finite-enumeration proof of the unbounded minimality statement.

The checker was independently written from the displayed mathematical formulas;
the source verifier was read but was not copied or executed. The report does not
redistribute the source TeX. SymPy is used as an external dependency.
