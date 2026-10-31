# Bala's A111984 congruence modulo five

For the normalized integer formal cube root

$$A(x)^3 = \sum_{k\ge0}(2k+1)8^k x^{k(k+1)/2},\qquad A(0)=1,$$

Peter Bala's FORMULA conjecture in [OEIS A111984](https://oeis.org/A111984/internal) states that $5\mid a(5n+4)$ for every $n\ge0$.
The sequence is due to Paul D. Hanna (2005); Bala's conjecture is dated February 26, 2021. The source snapshot is revision 18.

Run `python3 docs/reports/bala-a111984-mod-five/check.py` from the repository root. The program uses only exact integers and standard Python libraries. Its output is in [check.out](check.out).

The check matches all 18 displayed terms and obtains the coefficients through degree 1104 by two independent recurrences. It recomputes the defining cube and checks the truncated Jacobi identity, the congruence used for integrality, the Frobenius product identity, and both zero residue classes modulo five. All 221 tested indices in each of classes 3 and 4 pass. The first coefficient in each of classes 0, 1 and 2 provides a negative control against claiming extra zero classes.

These finite calculations support the proposed unbounded proof; they do not establish that proof. No Lean proof is included.

The integrality ingredient is the postscript of N. Heninger, E. M. Rains and N. J. A. Sloane, [On the Integrality of n-th Roots of Generating Functions](https://arxiv.org/abs/math/0509316v4), J. Combin. Theory A 113 (2006), 1732–1745, DOI [10.1016/j.jcta.2006.03.018](https://doi.org/10.1016/j.jcta.2006.03.018).

The 18 source terms are reproduced with attribution from OEIS A111984 under [CC BY-SA 4.0](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement). This report and its independently written program are distributed under the same license.
