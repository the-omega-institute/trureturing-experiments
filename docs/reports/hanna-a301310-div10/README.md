# Hanna square-binomial geometric sums

`python3 check.py` uses only Python 3 standard-library integer arithmetic. The
stored output is `output.txt`; success ends with `ALL_CHECKS_OK`. This is finite
support for preregistration, not a universal proof or a Lean verification.

The source is Paul D. Hanna's OEIS entries:

| Entry | Exact coefficient formula | Published conjecture |
|---|---|---|
| [A301310](https://oeis.org/A301310/internal), revision 22 | sum over k >= 0 of 2^k binomial(k^2,n)/3^(k+1) | divisible by 10 for n > 0 |
| [A301311](https://oeis.org/A301311/internal), revision 14 | sum over k >= 0 of 2^k binomial(k^2+n-1,n)/3^(k+1) | divisible by 10 for n > 0 |
| [A301432](https://oeis.org/A301432/internal), revision 7 | sum over k >= 0 of 3^k binomial(k^2,n)/4^(k+1) | divisible by 21 for n > 0 |

The three b-files are Paul D. Hanna's tables for n = 0..100, downloaded from
`https://oeis.org/A<number>/b<number>.txt`. Attribution and reuse are subject to
the [OEIS End-User License Agreement](https://oeis.org/LICENSE).

For each positive index n <= 200, put d = 2n and take the integer values of
the source polynomial P(k). One exact evaluation is the Newton sum
`sum_j Delta^j P(0) m^j`, where m = 2 or 3. A separate convolution computes the
coefficients of `(1-z)^(d+1) sum_k P(k) z^k`. The program compares their
homogeneous evaluations, checks two vanishing tail coefficients, coefficient
reciprocity, division by z(z+1) over the integers and the extracted factor
m(2m+1). Every source-table value through n = 100 is reproduced. The program
uses the standard finite polynomial/geometric-series identity; it does not
approximate an infinite sum by a finite partial sum.

Results: 303 exact source-table matches; 600 positive-index divisibility and
factorization checks. The output includes initial values, initial quotient
polynomials, digit counts and hashes of the computed sequences. The cases
n = 0 are checked for fidelity and excluded from the divisibility claims.

The neighboring entries A301312, A301466 and A301468 contain no published
COMMENTS divisibility conjecture and are not preregistration members here.
