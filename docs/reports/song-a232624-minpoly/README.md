# Song's OEIS A232624 coefficient gap

Primary source: [OEIS A232624](https://oeis.org/A232624), revision 109,
Jianing Song's July 11, 2019 conjecture. The coefficient-array examples use
the monic minimal polynomial of `2*cos(2*pi/n)` over the rationals.

Run `python3 check.py`. Python 3 standard library only; all polynomial
arithmetic uses integers. `output.txt` contains the result.

The program reproduces the entry's coefficient-table rows 1 through 20 and
its full half-difference example at n=15. For every odd n from 3 through 451,
it checks nonzeroness, the exact degree and leading coefficient, reflection,
and every stated coefficient consequence. It builds cyclotomic polynomials
by exact division and extracts real polynomials by triangular reciprocal
substitution. An independent Chebyshev divisor construction checks both
polynomials in every comparison.

These finite checks support source fidelity and detect small counterexamples.
They do not prove the unbounded conjecture and are not a Lean probe.
