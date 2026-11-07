# Zhai's small-function Conjecture 1

Source: Wenguang Zhai, *On Some Sums Involving Small Arithmetic Functions*,
Acta Mathematica Sinica (English Series) 40 (2024), 2497–2518,
DOI 10.1007/s10114-024-2129-y; arXiv:2607.17474v1, Section 1.2.

Run `python3 check.py`. The program uses Python's standard library and exact
integers and rational numbers. Its output is in `check.stdout`.

Fidelity checks reproduce the two exact divisor-convolution identities immediately
before Corollary 2 for all integers 1 through 1024, and check equation (1.1) for the
constant function on a grid of rational arguments. Independent direct enumeration
checks the quotient-count identity for X through 128 and the hole count for every
even T through 256.

For the proposed counterexample, T starts at 1024 and is repeatedly squared.
The fixed function is f(m)=1+h(m), where h indicates the integers in [T,2T) with
floor(T²/m)=floor(T²/(m+1)) for some such T. The program enumerates the first two
blocks, obtains exactly 512 and 524288 holes, and checks the interference and
endpoint conditions using integer arithmetic.

The infinite constant C_h is retained symbolically. Exact integer division onto
a rational lattice gives outward enclosures for the two enumerated block masses;
the entire remaining infinite tail is bounded by the telescoping majorant 1/T₂.
The full errors are exactly S_h(X)-X C_h, with S_h(X) computed as an integer.
The printed rational enclosures certify E_f(X)<-sqrt(X)/16 at both checked scales.
These are rigorous enclosures, not exact finite evaluations of the infinite series.

The all-scale conclusion and the exclusion of every exponent below 1/2 require
the preregistration's general argument. This finite check does not prove them,
and does not address or refute the source's mean-square theorem.

The program is original code for this report; the repository's Apache-2.0 licence applies.
