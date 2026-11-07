# Ehrenborg odd-D pizza coefficient checks

Source: Richard Ehrenborg, *Conjectures for cutting pizza with Coxeter arrangements*, arXiv:2410.13593v1, Conjecture 5.2 and Table 2. The source's version is fixed; these checks do not certify the universal conjecture.

Run with Python 3 and SymPy:

    python3 check.py > check.out

The arithmetic uses integers, `fractions.Fraction`, and exact SymPy polynomials. The program:

- reproduces both displayed Table 2 coefficients for each of D3, D5, D7;
- compares direct rotated odd-power coefficients with the proposed Gamma-moment identity for every p,q from 0 through 12;
- verifies negativity and equality of the two coefficient formulas for all partitions of weight at most 12, 6, 3, 2 in dimensions 3, 5, 7, 9 respectively;
- verifies the corrected coefficient projection for the geometric consequence at n=3, d=6, where the unprojected specialization differs by 7(a^6-b^6)/720;
- checks the reciprocal-distance Pfaffian sign on every subset of {1,...,8} of sizes 2, 4, 6;
- cancels the determinant times the reciprocal-distance Pfaffian to a polynomial before applying product moments in four exact cases.

The direct calculation expands the source's two rotated odd powers and factors its multinomial sum over disjoint pairs. It does not use the proposed moment formula. Thus equality of the two implementations checks normalization, matching signs, and factorial factors independently. The factor relating a leading coefficient to Z_d/J is 2^(n(n-1)/2).

`check.out` is the complete output; success requires exit 0 and the final line `ALL_CHECKS_OK`. The experiment is finite evidence only. A settlement needs a general proof, including the coefficient identity, strict Pfaffian positivity, and cancellation before integration.
