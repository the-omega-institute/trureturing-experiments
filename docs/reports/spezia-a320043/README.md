# Exact checks for Spezia's A320043 square-value conjecture

Source: [OEIS A320043](https://oeis.org/A320043/internal), revision 43,
sole COMMENTS conjecture and first FORMULA; the linked table is
[OEIS A322550](https://oeis.org/A322550/internal).

For n ≥ 1, the source defines

    A(n) = sum(k=1..n, (n+1-k)^2*k/gcd(n+1-k,k)^3).

Run with Python 3 and SymPy (checked with Python 3.14.8 and SymPy 1.14.0):

    python3 docs/reports/spezia-a320043/check.py

The saved output is `check.out`. All evidence uses exact integers, rational
arithmetic, or symbolic polynomial expansion. The program:

- reproduces all 42 published DATA terms and A(96) = 7376656;
- compares integer and rational evaluation for every 1 ≤ n ≤ 256;
- verifies the polynomial telescope and the prime-row identity directly for
  all 55 primes p ≤ 257;
- checks the Pell doubling chain, the complete Lucas certificate for
  p = 708158977, and independent trial division through 26611;
- verifies the exact proposed witness n = 708158976 with square root
  144767444036350576 and square 20957612852817895952216852375531776.

The large row is inferred from the prime-row identity, whose difference and
endpoints are checked symbolically; its 708158976 summands are not evaluated.
The accompanying mathematical argument must justify the gcd reduction and
summation. These checks are experimental evidence for a proposed refutation;
no Lean proof or kernel-checked settlement is supplied here. Infinitude of
prime Pell solutions and completeness of the remaining square values are not
asserted.
