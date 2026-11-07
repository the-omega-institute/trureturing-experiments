# Voloshyn boundary factorization: exact bounded checks

Source: S. Voloshyn, *Differential Recursions, Projective Discriminants, and Algebraic Generating Functions*, arXiv:2609.31189v1, Conjecture 8.1, equations (4)–(10), (95)–(97), and the displayed p=4 and specialized p=6 boundary formulas. Source URL: https://arxiv.org/abs/2609.31189v1.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The stored stdout is `check.out`; success requires exit 0 and final line `ALL_CHECKS_OK`.

For q=2,...,6 the program:

- reproduces the source's R1 and R2 rational functions directly from its differential recurrence;
- compares recurrence coefficients with the proposed closed coefficient formula for k=0,...,3 and x-degrees 1,...,8;
- solves the proposed polynomial branch at x=1,2,3 through degree four in s, using exact rational arithmetic, and compares the transformed f coefficients with the recurrence;
- checks the rational primitive identity, its boundary evaluations, and the two reduced rational-map degrees;
- computes the complete boundary norm by an exact resultant with symbolic x and compares it with the source's product over roots of unity, also evaluated by an exact resultant.

The program additionally reproduces the source's displayed p=4 boundary polynomial and the p=6 boundary at x=1, and checks the separately normalized Narayana endpoint. It uses integers, Python fractions, and exact SymPy algebra throughout. The output supports the proposed ordinary argument; these bounded checks neither prove the all-order conjecture nor certify minimality of the norm, and they are not Lean verification.

The program is independently written for this report. The source supplies the recurrence and conjectured boundary formulas; the proposed auxiliary rational maps are

U(Z)=(1-Z)^q(1-y^q/Z^q),

V(Z)=(1-Z)^(q-1)(1+(q-1)Z-q y^q/Z^(q-1)),

f(Z)=(-1)^(q+1)(V(Z)-(1-y)^q)/q^3,

with y=-x and s=(-1)^q q^2 H. No source text or external implementation is vendored.
