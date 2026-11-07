# Pointwise periodic polynomial entropy: exact auxiliary checks

Run `python3 check.py` using Python 3. The program needs no third-party packages.
All arithmetic uses integers and `fractions.Fraction`. The saved output is
`check.out`; a successful run ends with `ALL_CHECKS_OK` and exits zero.

The source is Maša Đorić, Jelena Katić and Milan Perić, *A flexibility result for
polynomial entropy of pointwise periodic homeomorphisms*,
[arXiv:2606.28283v1](https://arxiv.org/abs/2606.28283v1), Sections 1–3 and the
unnumbered closing Question in Section 4.3.

The fidelity check reproduces the source's periodic-homeomorphism zero-growth
result on a five-point, non-isometric system with cycles of lengths two and three.
An exact period-six certificate implies that its Bowen metrics stabilize at
horizon six. Its separated cardinalities therefore stabilize. The program also
checks all four source covering/separation quantities in 90 finite instances.
It uses the standard separated-set convention: the inequality is required for
distinct points. The source's displayed prose omits this qualification.

The proposed affirmative settlement uses the planar family

\[
a_{n,k}=2+\cos(2\pi k/n),\quad
\beta_{n,k}=2^{-n}(1+k/n),\quad
v_{n,k}=a_{n,k}(1,\beta_{n,k}),\qquad n\ge2,\quad0\le k<n.
\]

Join each endpoint to the origin, add the fixed interval \([0,3]\times\{0\}\),
and cyclically permute each star, preserving its radial fraction. The target
ordinary argument gives covering dimension one and polynomial entropy one half.
The program checks the following auxiliary facts:

- Exact disjoint positive slope blocks and tail-height bounds.
- Modular observation-time certificates separating endpoints of different periods.
  The reduced phase belongs to \([1/3,1/2]\); the ordinary analytic proof uses
  the nonpositivity of cosine there. Eighteen actual Euclidean squared-distance
  checks use only rational cosine values, for three nondegenerate choices of
  \(c>r>0\) in the profile \(c+r\cos\).
- The strict observation horizon, including horizons just below, at and above
  squares.
- Rational overestimates for the parameter-cell diameters and the number of
  reciprocal-frequency cells. The analytic inequality \(2\pi<8\) supplies the
  Lipschitz constant; the program does not infer that inequality from samples.

These checks are finite auxiliary evidence. They do not establish compactness,
covering dimension, continuity, the infinite-family entropy limits, literature
novelty, or a Lean proof. Those remain obligations of the general argument and
subsequent formalization. The positive settlement's main content is that general
argument; these finite checks are not its admission basis.

Closely related cyclic finite-star geometry appears in Daghar–Naghmouchi,
[arXiv:2109.04246](https://arxiv.org/abs/2109.04246), Example 4.9. Its hyperspace
entropy calculation was corrected in Naghmouchi,
[arXiv:2512.08836v3](https://arxiv.org/abs/2512.08836v3), Remark 2.4(ii).
Neither of those passages computes the base-map polynomial entropy used here.
The new candidate contribution is the quantitative separation and covering
argument, rather than the finite-star motif.
