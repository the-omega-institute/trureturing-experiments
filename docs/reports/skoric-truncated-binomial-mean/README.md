# Exact checks for Škorić's truncated binomial mean

Source: Boris Škorić, *Quantum digital signatures with smaller public keys*,
arXiv:2012.15493v2, Conjecture 1, equation (35); Quantum Information and
Computation 21 (2021), 955–973, DOI 10.26421/QIC21.11-12-4.

For positive integer caps, the source defines

\[
f(n,r)=\frac1r\frac{\sum_{k=0}^r k\binom nk}{\sum_{k=0}^r\binom nk},
\qquad r\le n/2,
\]

and conjectures decrease in the cap. Run `python3 check.py`; the expected output
is in `output.txt`. The program uses integers and `fractions.Fraction` throughout.
It requires Python 3 and no third-party package.

The fidelity check reproduces the exact equality in source equation (45),
\(f(2z,z)=(1+\binom{2z}{z}/2^{2z})^{-1}\), for every integer \(1\le z\le100\).

The settlement checks cover every pair \(1\le r<s\le\lfloor n/2\rfloor\)
for \(2\le n\le200\), and every adjacent pair with \(1\le r<n\le200\).
They also check the proposed proof's algebraic identity and convolution-block
inequalities on 786 finite positive ultra-log-concave rows. These rows are
constructed by taking the ratios \((k+1)a_{k+1}/a_k\) nonincreasing from
\(\{1/3,1/2,1,3,8\}\), with 2 through 7 ratios.

Separate extension checks cover 305 rational upper arguments in \([2,40]\),
using falling-product binomial coefficients and positive integer caps at most
half the upper argument, and 116 interior Poisson-binomial laws with 2 through
7 coordinates in \(\{1/5,1/2,4/5\}\), up to coordinate permutation.
These extensions are not additions to the registered source claim.

All comparisons are exact, including the boundary cap 1. A failed assertion
exits nonzero. Success prints `ALL_CHECKS_OK` as its final line. These finite
checks do not prove the conjecture or certify the proposed universal argument.
No Lean proof is included.
