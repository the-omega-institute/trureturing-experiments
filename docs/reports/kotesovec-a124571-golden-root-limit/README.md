# OEIS A124571: exact checks for the golden root limit

`python3 check.py` uses only Python 3 standard-library exact integers and rational numbers. Its expected output is `check.out`, ending in `ALL_CHECKS_OK`.

The source is [OEIS A124571, revision 20](https://oeis.org/A124571/internal), whose COMMENTS contain Hanna's root-limit question and Kotesovec's GoldenRatio conjecture. The 38 embedded terms reproduce the complete displayed sequence prefix. Source terms and definitions are attributed to the OEIS contributors; OEIS material is available under [CC BY-SA 4.0](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement). The independently written check program is released under CC0-1.0.

The program computes coefficients through degree 512 using Hanna's FORMULA recurrence. It independently checks the original b-recurrence, the normalization, monicity and full interval support, and the critical identity in the exact quadratic ring where rho squared equals one minus rho. It builds balanced trees and verifies their binary-digit node counts and exact index-gap identities. Polynomial coefficient comparisons check the tree lower bound for N=1..3 and m=1..24. Rational brackets for four finite roots are certified by integer powers.

The checks support a proposed general analytic and combinatorial proof. A finite experiment cannot verify existence or the value of an infinite limit, and the report does not claim otherwise. No Lean proof or refutation is included.
