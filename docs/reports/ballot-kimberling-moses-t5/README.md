# Ballot–Kimberling–Moses Table 5: the double-root row

The target is the polynomial-form tree modulo `q(X) = X^2 - 2X + 1` in Christian Ballot, Clark Kimberling and Peter J. C. Moses, *Linear Recurrences Originating From Polynomial Trees*, The Fibonacci Quarterly 55(5) (2017), 15–27, Section 5, Table 5, p. 23. DOI: [10.1080/00150517.2017.12427731](https://doi.org/10.1080/00150517.2017.12427731). The [published PDF](https://www.fq.math.ca/Papers1/55-5/BallotKimberlingMoses.pdf) supplies the reference values.

Run `python3 docs/reports/ballot-kimberling-moses-t5/check.py` from the experiments root. Python's standard library is sufficient. The program uses integers throughout; `check.out` is its complete output.

The fidelity checks reproduce Theorem 1.1(1), the sizes of all generations 0 through 14 of the universal polynomial tree, and verify Theorem 1.1(4)'s degree-plus-coefficient-sum generation formula for its 16,383 nonzero nodes. They also reproduce every displayed value of the selected Table 5 row, G_0 through G_10.

Two separate breadth-first searches enumerate the quotient tree through depth 100. One computes polynomial remainders by integer long division; the other uses the proposed coordinates `(a,b) = (p(1),p'(1))` and successors `(a+1,b)`, `(a,b+a)`. Their complete shells agree. For all 166,751 visited residues, the first-generation formula is zero at `(0,0)` and `a + ceil(b/a)` elsewhere. The counts equal `1 + binomial(n,2)` in all 101 tested generations, and the proposed recurrence holds at every tested index 3 through 100. The output ends with `ALL_CHECKS_OK`.

These finite checks supply fidelity and counterexample screening. They do not establish the formula or recurrence at unbounded depth. The planned positive settlement requires a general proof of the exact distance and shell count. No Lean proof is included.
