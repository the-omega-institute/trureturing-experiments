# Spezia's A322844 serpentine matrix conjectures

`check.py` checks the two conjecture lines in [OEIS A322844](https://oeis.org/A322844), revision 68, with the matrix definition from [A317614](https://oeis.org/A317614), revision 194. The entries are the integers 1 through n² written alternately forward and backward along successive rows.

Run `python3 check.py` with Python 3 and SymPy. The recorded output uses Python 3.14.8 and SymPy 1.14.0. All evidence uses integer or exact rational arithmetic.

Fidelity checks reproduce M(2), M(3), the first ten A317614 trace values, twenty A322844 values, and the source's proved determinant vanishing for n > 2. For 0 ≤ n ≤ 40, the program independently computes the full characteristic polynomial, all principal minors of sizes two and three, exact matrix rank, and the proposed rational factorization and four finite sums. It also evaluates every principal minor of size at least three for n ≤ 10.

The second exterior trace equals the negative of the defining quartic; the first four signed values are −5, −6, −68, −50. Every tested higher exterior trace vanishes. The output ends in `ALL_CHECKS_OK` and the check exits 0. Finite checks support preregistration; they do not prove either universally quantified conjecture.
