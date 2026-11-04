# Han–Molnár trace-centrality checks

The target is the general C*-algebra question immediately after Lemma 11 in
Jiayi Han and Lajos Molnár, *Notes on trace characterizations via inequalities
concerning different versions of the geometric mean and relative entropies*,
J. Math. Anal. Appl. 567 (2027), 131130, DOI
[10.1016/j.jmaa.2026.131130](https://doi.org/10.1016/j.jmaa.2026.131130).
The institutional published-version text is available at
https://publicatio.bibl.u-szeged.hu/41664/1/1-s2.0-S0022247X26007420-main.pdf.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. Its output is in
`check.out`. All arithmetic uses integers, `fractions.Fraction`, or SymPy exact
algebra. No floating-point approximation supplies evidence.

The source-fidelity check reproduces the rank-one geometric-mean formula in
the proof of Lemma 11 (printed pp.20–21) in 720 exact cases. It evaluates the
matrix fractional power independently using SymPy's exact spectral algebra.

The settlement checks test the corner-inverse identity `X C X = X` and its
projection consequence for all 300 selected two-dimensional cases, including
complex projections, and eight rank-two cases in dimension three. The
geometric mean equals X at three rational parameters in every two-dimensional
case and at parameter 1/2 in the three-dimensional cases. The selected
Schur-resolvent gaps are checked for positivity and equality exactly.

An additional rational interval certificate verifies a strict violation of
the original inequality at invertible X: parameter 1/2, B = diag(1,1/4), and
X = (1/5)[[1,1],[1,1]] + (1/10000)I. With H equal to right minus left after
cancelling B-X, the program proves -1/10 < H < -1/20. Logarithm bounds use a
positive atanh series with an explicit rational tail bound; square-root
bounds use integer square roots.

These checks support the proposed general proof mechanism. They neither
prove the general C*-algebra implication nor establish novelty or literature
completeness. No Lean proof was written. The paper's text and PDF are not
redistributed here.
