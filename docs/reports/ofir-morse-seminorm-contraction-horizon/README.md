# Ofir–Morse seminorm contraction horizon

The independent Python 3 program uses only integers and `fractions.Fraction`.
Run `python3 check.py` from this directory; the stored output is `check.out`.
A successful run exits 0 and ends with `ALL_CHECKS_OK`.

Source: Ron Ofir and A. Stephen Morse, *On Determining the Convergence Rate of
an Infinite Product of Stochastic Matrices*, arXiv:2607.03623v1, submitted
3 July 2026, Sections I–II, IV Corollary 2, and V's final paragraph.
https://arxiv.org/abs/2607.03623v1

The fidelity check reproduces Section II Proposition 1's equality between the
Dobrushin coefficient and the induced quotient infinity seminorm on three
specific stochastic matrices. It evaluates the quotient objective at all eight
cube corners independently of the row-difference formula, including the
non-scrambling boundary matrix P. A further instance of Corollary 2 certifies
eventual contraction of A_1 using the rational Frobenius-square upper bound at
power 512.

The candidate family uses J=11^T/3, H=I-J,
P=[[1,0,0],[0,1,0],[1,0,0]], T=HP, S=(29/30)P+(1/30)J,
and A_K=(1-1/(2K))I+S/(2K). The fixed seminorm under consideration is
nu(M)=||HM||_2 on the equal-row-sum subspace; the norm is the Euclidean operator
norm. This definition is independent of K. The program checks stochasticity,
strict positivity, stationarity of pi=(59,30,1)/90, projection identities,
the formula for H A_K^j, and the bilinear certificate at x=(2,-3,1).
Since x^T x=14, Cauchy–Schwarz makes x^T H A_K^j x>14 a certificate
for nu(A_K^j)>1. No numerical singular-value calculation is used.

The tested horizons are 1 through 12, 20, 50 and 100, with every positive power
through each horizon checked (248 matrix powers). Exact scalar bounds are also
checked for j=0 through 100 on the listed rational grids. These finite checks
validate arithmetic and supply regression data. They do not establish the
universal quantifier over K, the seminorm interface, or the limit proof; those
remain obligations of the proposed mathematical settlement.

The program is original code under this repository's license. It copies no
third-party implementation. The cited mathematical source supplies the
comparison statements; no source TeX is redistributed here.
