# Kotesovec weighted distinct-partition logarithmic asymptotics

Run `python3 check.py` using Python 3 and its standard library. `output.txt` records a successful run ending in `ALL_CHECKS_OK`. The program uses integers and `fractions.Fraction`; it uses no floating point.

For each r in {1,2,3,4,5}, the coefficient sequence is defined by the formal product `product_{k>=1}(1+k^r*x^k)`. The conjectured logarithmic ratio is `log(a_r(n))/(r*sqrt(n/2)*(log(2*n)-2)) -> 1`.

The checker reproduces all displayed initial terms of OEIS A022629, A092484, A265840, A265841 and A265842, respectively 43, 31, 29, 26 and 23 terms, and the two worked examples at n=6. Recursive distinct-partition enumeration independently checks all five weighted-subset identities for n=0..45. Dynamic programming checks the triangular witness, the integer AM-GM inequality and the count-times-maximum bound for n=1..180. Rational logarithm and square-root enclosures certify the counting bound and the final proposed logarithmic sandwiches for every n=2..180 and r=1..5. The logarithm enclosures use the positive atanh series with a geometric tail bound and outward rational rounding; square roots use integer `isqrt`.

These are finite checks supporting a proposed paper outline. They do not prove an infinite limit, settle an OEIS conjecture or certify a Lean proof.

Sources: [A022629](https://oeis.org/A022629), revision 54; [A092484](https://oeis.org/A092484), revision 31; [A265840](https://oeis.org/A265840), [A265841](https://oeis.org/A265841), [A265842](https://oeis.org/A265842), revision 11 each. The conjectures are attributed to Vaclav Kotesovec. OEIS sequence data and quotations are subject to the [OEIS End-User License Agreement](https://oeis.org/LICENSE). The initial-term lists in `check.py` are attributed to those records. The check program is an independent implementation.
