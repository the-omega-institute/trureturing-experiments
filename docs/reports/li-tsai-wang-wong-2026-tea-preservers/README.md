# Exact checks for complex TEA preservers

Source: Chi-Kwong Li, Ming-Cheng Tsai, Ya-Shu Wang and Ngai-Ching Wong, “Linear maps of C*-algebras preserving norm parallelism,” [arXiv:2611.06187v1](https://arxiv.org/abs/2611.06187v1), Problem 1.1, second sentence.

Run with Python 3 using `python3 check.py`. The program uses only the standard library, Gaussian integers represented by integer pairs, and `fractions.Fraction`. Its expected last line is `ALL_CHECKS_OK`; the complete output is in [check.out](check.out).

The fidelity check reproduces Lemma 2.1(a.2) at beta = 1/2. On 81 Gaussian input vectors it verifies all 6,561 parallel-pair equivalences for the matrix [[2,1],[1,2]], whose scale is irrelevant to both relations. It also exhibits an original TEA pair whose image squared norms are (5,5,16), so the images fail triangle equality. This confirms the source's distinction between complex parallelism and TEA preservation.

The settlement screen tests every two-row map with nine Gaussian integer choices per coefficient (6,561 maps), every one-row map with coefficients in {0,1,-1,i,-i} (25 maps), and every three-row map with these five choices (15,625 maps). Each map is tested on all 233 TEA pairs from the 25 input vectors with coordinates in {0,1,-1,i,-i}. Every surviving map is checked separately to have coordinate-supported rows, equal maximum column scales, and the predicted scalar sup norm on all tested inputs and their sums. There are respectively 65, 1 and 481 survivors. The one-row survivor is zero. These maps also test seminorms into arbitrary finite sup-norm codomains of the indicated dimensions.

For nonnegative squared norms A, B, C, equality sqrt(C) = sqrt(A) + sqrt(B) is tested exactly by C-A-B >= 0 and (C-A-B)^2 = 4AB. No square-root approximation or floating point tolerance enters any decision.

A separate rational calculation checks the phase perturbation in the proposed two-coordinate seminorm lemma: alpha = 1, beta = 1/4, s = -16 and t = 1/16 give squared endpoint values 256/257, below the squared midpoint value 1. A finite-support c_0 model checks the rational local cutoffs and all nine corners of the nonunital boundedness argument.

These are finite fidelity and consistency checks. They prove neither the assertion for arbitrary complex C*-algebras nor the general two-coordinate seminorm lemma. Pure-state norm attainment, irreducible GNS representations, Kadison transitivity and the hereditary-subalgebra limit argument are ordinary-mathematical ingredients of the preregistration outline; they are not certified by this program.
