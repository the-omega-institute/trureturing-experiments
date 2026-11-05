# Square chessboard cover ideals: exact checks

Source: Arka Ghosh and S. Selvaraja, *Powers of Alexander Duals of Stanley-Reisner Ideals of Chessboard Complexes*, [arXiv:2609.37219v1](https://arxiv.org/abs/2609.37219v1), Theorem 4.7, Corollary 4.9 and the first question of Question 4.10.

Run `python3 check.py` with Python 3. The program uses only the standard library, integers and `fractions.Fraction`. The stored result is `check.output.txt`; its final line is `ALL_CHECKS_OK`.

The program forms the minimal ordinary generators from permutation complements and the minimal symbolic generators from the edge exponent inequalities. For boards of sizes 1 and 2 it computes graded Betti numbers by tensoring the Taylor resolution with the residue field and performing exact row reduction, separately over the rationals and the fields with 2 and 3 elements. It reproduces the source's square-board regularity formula at size 2 and exponents 2 through 6, and its comparison corollary at exponents 1 through 6. The requested comparison holds in these computations at sizes 1 and 2 and exponents 3 through 6.

For sizes 2 through 4 and exponents 1 through 4, it checks that the proposed lcm has exactly the two claimed ordinary-power generators below it and that every immediate proper divisor has at most one. For sizes 2 and 3 and exponents 2 through 4 it checks the colon identity for every independent set and all exponent vectors in the stated window: 1,300,239 cases. This window tests membership thresholds after the colon; larger individual exponents saturate the same inequalities. Monomials involving a variable in the added-variable ideal belong to both sides. It also checks the induction arithmetic for sizes and exponents at most 30.

The triangle example confirms that symbolic and ordinary squares of a clique cover ideal can differ. That assertion occurs in the source's separate rectangular size-1 argument. The proposed square-only induction does not use it.

These computations support a general argument; they do not prove Question 4.10, field independence beyond the tested fields, or any Lean homological bridge. The broader equality classification in Question 4.10 is outside this report.
