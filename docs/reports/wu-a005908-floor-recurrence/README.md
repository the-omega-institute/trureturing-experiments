# Wu's A005908 floor recurrence

Source: [OEIS A005908](https://oeis.org/A005908/internal), revision 27
(15 January 2022); Chai Wah Wu's two FORMULA conjectures (13 June 2020).
The source's Mathematica generator fixes offset zero and initial values 0, 1.

Run `python3 docs/reports/wu-a005908-floor-recurrence/check.py` from this
repository with Python 3 and SymPy. The program uses integer arithmetic and
exact symbolic algebra, requires no network, and prints `ALL_CHECKS_OK` only
after every check succeeds. Its output is stored in `output.txt`.

The fidelity check reproduces all 30 displayed OEIS terms. Further checks
compare the original floor recursion, the integer companion and the rational
generating-function coefficients through degree 3999. The symbolic checks
verify the factorization, companion relation, centered error recurrence,
contracting quadratic energy and its strict positive margin. The finite
comparisons support the proposed argument for every index; these computations
do not establish either conjecture for all indices and contain no Lean proof.
