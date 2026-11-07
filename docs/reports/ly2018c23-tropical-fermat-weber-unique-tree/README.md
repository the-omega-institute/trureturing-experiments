# Exact checks for Lin–Yoshida Conjecture 23

Run with Python 3; no external packages are needed:

```sh
python3 check.py
```

The program reproduces the three-point example in Proposition 10 of Bo Lin and Ruriko Yoshida, *Tropical Fermat-Weber points*, arXiv:1604.04674v3, DOI 10.1137/16M1071122. It then checks a five-leaf, four-observation nonstar witness, 110 generated positive binary trees with seed 337, and all 1,028 partitions of 3 through 7 coordinates into at least three nonempty classes.

Each case has a balanced affine lower-bound certificate for global optimality. The exact equality constraints form a strongly connected directed graph, certifying uniqueness throughout the ambient quotient by constant vectors. Strict deletion competitors certify essentiality. For tree cases every ordered ultrametric inequality, rooted triplet, and positive common-root-height normalization is checked. The two equivalent metric formulas are compared exactly.

`output.txt` is the output of `check.py`; success ends in `ALL_CHECKS_OK`. All evidence uses integers and `fractions.Fraction`. The finite checks do not establish the universal family theorem or a Lean kernel result.

The check program is original repository experiment code under the experiments repository's licence. The paper supplies the mathematical definitions and the fidelity example; no third-party program source is copied.
