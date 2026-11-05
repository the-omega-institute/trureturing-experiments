# Hanna A266481 prime-multiple congruence: exact finite checks

The source is Paul D. Hanna's unnumbered conjecture in [OEIS A266481](https://oeis.org/A266481), revision 35, with connected logarithmic coefficients from [A266526](https://oeis.org/A266526), revision 9, and the polynomial triangle [A266521](https://oeis.org/A266521), revision 13. The conjecture says a(p*n) is one modulo p for every prime p and natural n.

Run `python3 docs/reports/hanna-a266481-prime-congruence/check.py` from the repository root. Python 3's standard library suffices. Integers and fractions provide exact arithmetic; no numerical approximation is used. The saved output is `check.out`, whose final line is `ALL_CHECKS_OK`.

The program reproduces all 23 displayed A266481 terms, all 20 displayed A266526 terms, and the first six A266521 rows. It computes the defining formal-root limit directly through degree six using Laurent polynomials in t=1/N. Independently, Prüfer tree enumeration with orientations and injective assignments to two distinguished outgoing slots agrees with the connected coefficients through six vertices. Enumeration of all slot forests agrees with the sequence through five vertices. The independent block-rotation actions preserve these forests and fix only the empty forest for (p, blocks)=(2,1), (2,2), (3,1), (5,1). Reciprocal arrows on two vertices are rejected as a negative control. All 53 prime/index pairs with p<=30 and p*n<=30 satisfy the conjecture.

A266521's example for exp(A(x,y)) has a first-term typo: it prints (1+y)*x, while its defining NAME and its displayed first connected row give (1+y)^2*x. The checks use the defining NAME and verify the displayed connected rows. A266481's independent formal limit agrees.

These finite checks support source fidelity and the proposed forest interpretation. They do not prove the conjecture or the all-degree counting identity. The logarithmic weight w_4=376/3 also records why the integral-weight hypothesis of Ahaan Kallat's arXiv:2610.15849v1 cannot be applied to this sequence as stated.

The numerical source data and definitions are credited to the OEIS contributors and remain subject to the [OEIS license](https://oeis.org/LICENSE). The independently written verification program follows this repository's license.
