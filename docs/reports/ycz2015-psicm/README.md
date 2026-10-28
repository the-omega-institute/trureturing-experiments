# Yang–Chu–Zhang psi complete-monotonicity certificate

Source: Zhen-Hang Yang, Yu-Ming Chu and Xiao-Hui Zhang, *Necessary and sufficient conditions for functions involving the psi function to be completely monotonic*, Journal of Inequalities and Applications 2015:157, DOI [10.1186/s13660-015-0674-8](https://doi.org/10.1186/s13660-015-0674-8), Conjecture 1(ii). Earlier source: Yang, arXiv:1409.6411v1.

Run with Python 3 and SymPy (`python3 check.py`). All constants use integers, `fractions.Fraction` or exact SymPy rationals. `output.txt` is the program output. No floating point is used as evidence.

The fidelity check reproduces the exact rational identity in Remark 4 following Corollary 4 and verifies its negative sign for every integer n ≥ 1 by positive coefficients after shifting n = m+1. It also cross-checks n = 1,…,32.

The settlement check rebuilds

\[
R(x)=\psi(x+1/2)-\log x,\quad
G_a(x)=24(x^4+ax^2+7a/40-31/336)R(x)-(x^2-7/40+a).
\]

For A = 155/294, ε = 1/10000 and r = 1/64, every a ∈ [A−ε,A] is positive. Put u = r+1/2 = 33/64 and c = 7(a−A)/40. The analytic inputs are the classical convergent digamma series (DLMF 5.7.6), its derivatives on compact subsets of u > 0, γ < 1, the atanh series for logarithms and the integral comparison for a decreasing positive function. These facts are mathematical inputs, not established by executing Python.

The independent bounds are

\[
R(r)>-1-u^{-1}+\sum_{k=1}^3(k^{-1}-(k+u)^{-1})+6\cdot56/81
=394328567/231912450>5/3,
\]

\[
R'(r)>-64,\qquad
-\psi''(u)<2[u^{-3}+(u+1)^{-3}+1/(2(u+1)^2)]
=511623041024/32798729601<16<18.
\]

Here log 2 > 2(1/3+(1/3)^3/3) = 56/81. The series tails discarded in the first inequality are positive; the last tail is bounded by the integral from 1 to infinity. Consequently R''(r) > 4078. The program checks the differentiated product identity symbolically, all coefficient signs on the entire parameter interval, and the exact uniform estimate

\[
G_a''(r)>249085992523/481689600000>1/2.
\]

Thus the second derivative of −G_a is negative at an interior domain point throughout that interval, violating the n = 2 condition of the conjectured sufficient direction. This is an analytic certificate checked with exact algebra, not a Lean proof. It does not determine the corrected sharp cutoff or Conjecture 1(i).

Classical inputs: [DLMF §5.7](https://dlmf.nist.gov/5.7) and [§5.15](https://dlmf.nist.gov/5.15). The source article is open access under CC BY 4.0; this directory contains independently written code and derived exact data, with no copied article file.
