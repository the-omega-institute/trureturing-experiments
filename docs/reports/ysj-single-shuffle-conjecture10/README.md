# Single-shuffle Conjecture 10: exact expectation-of-norm checks

This report checks the proposed refutation of equations (25) and (26) in
Appendix B, Conjecture 10, of Chulhee Yun, Suvrit Sra and Ali Jadbabaie,
*Can Single-Shuffle SGD be Better than Reshuffling SGD and GD?*,
[arXiv:2103.07079v1](https://arxiv.org/abs/2103.07079v1).
It makes no Lean settlement claim.

Run with Python 3 and SymPy:

```sh
python3 docs/reports/ysj-single-shuffle-conjecture10/check.py
```

The checked output is in `check-output.txt`. Every assertion uses exact
symbolic algebra, integers or rational numbers. No floating-point sign test
or numerical approximation supplies evidence.

The program reproduces the Appendix A identity for the mean of permutation
product powers for K = 1 through 8. It then checks the six one-epoch products,
the 36 independent double products, their common determinants, the mean
squared Frobenius norms and the polynomial gap for the parameter q. A
nonnegative decomposition certifies the needed sign on 1/2 <= q <= 1.
The singular-value transfer identities and concavity derivatives are checked
symbolically on their proper interior domains. At q = 3/4, integer-square-root
rational enclosures directly certify that the repeated-product expectation
exceeds the independent-product expectation, for both the spectral norm and
its square, by more than 1/10^8.

The projectors originate in Appendix A of Yun–Sra–Jadbabaie. Their deformation
B_i = q I + (1-q) P_i and the displayed polynomial gap are also used in
Binghui Peng, *A Resolution of the SS–RS–GD Inequalities*,
[arXiv:2607.22620v1](https://arxiv.org/abs/2607.22620v1), Section 2 and the violation lemma (Lemma 6).
Peng's theorem concerns the norm of an expectation and uses a tensor square.
The proposed additional argument uses the fixed determinant of real 2 by 2
products and finite Jensen/Cauchy–Schwarz to compare expectations of norms.
The program is supporting exact evidence; the complete family argument and
its literature boundary belong to the preregistration, and formal proof is
subsequent work.

Validation environment: Python 3.14.8; SymPy 1.14.0.
