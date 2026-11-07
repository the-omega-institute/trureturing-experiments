# Odd orbit powers: exact factorization lengths

Source: J. A. Vega Coso, *The semiring of spectral duality skeletons:
cancellation, characters and factorization*,
[arXiv:2610.23684v2](https://arxiv.org/abs/2610.23684v2), Section 7,
item 4, first question; Theorem 4.18 and Remark 4.20.

Run with Python 3 and SymPy (checked with SymPy 1.14.0):

```sh
python3 check.py
```

The program uses integers, `fractions.Fraction`, exact polynomial factorization
and finite dictionaries. It reproduces Remark 4.20 for j = 1, 2, 3 and
exhaustively enumerates unordered atom factorizations for j = 1 through 5,
at t = 2 and t = 3. It also checks every factor in the proposed general
witness formulas at those parameters. `output.txt` contains the results.
Both parameters give length sets {2,...,2j+1}. At j = 5 the numbers of
unordered factorizations at lengths 2 through 11 are respectively
91, 39, 35, 12, 10, 3, 3, 1, 1, 1.

## Why the enumeration is exhaustive

Write a = t/2 > 0 and P = x^a + x^(-a). For any nonnegative-coefficient
factorization Q R = P^n, fix an exponent v in the support of R. The support
of Q plus v lies in the support of P^n, so differences of exponents of Q
are in 2a times the integers. The same holds for R. Translate each factor's
smallest exponent to zero. The translated polynomials have nonnegative
integer coefficients in z = x^(2a), with product (1+z)^n. Their leading
coefficients multiply to one. Unique factorization in the ordinary integer
polynomial ring therefore makes them (1+z)^k and (1+z)^(n-k).
Symmetry centers their exponents, recovering Q = P^k and R = P^(n-k).
This includes factors with arbitrary real exponents.

The profile P^k has constant term zero for odd k and binomial(k,k/2) for
even k. The program lists every allowed fixed-point count f between zero
and that term with the required parity. The only degree-zero pair is the
unit (1,1). It tests every proper binary split (k,f) = (i,g)(k-i,h), with
gh = f, to decide irreducibility. Since any divisor of such a pair is
also a divisor of the target, this list contains every relevant atom.
A recursion then enumerates every multiset of these atoms of total profile
degree n and fixed-point product zero, with no ordering duplicates.
SymPy independently checks the normalized polynomial factorization for
n = 1 through 11. Explicit convolution and all 2^k choices of orbit sites
independently verify the real-exponent profiles at both rational values of a.

The finite results support a preregistration with expected settlement
**Proved**. The proposed settlement uses a general argument, with one or two
f = 2 atoms and remaining single-orbit factors; the finite enumeration is
supporting evidence. No Lean proof or kernel-checked settlement is supplied.
