# Lorentz-entanglement annihilation and convexity

Run `python3 docs/reports/piana2506/check.py` from the repository root with Python 3.12 or later and SymPy. The expected transcript is in `check.out`.

The first test reproduces the strict subadditivity violation for the two matrices in Appendix B of La Piana and Müller-Hermes, *Annihilating and breaking Lorentz cone entanglement*. The remaining tests use the Appendix B criterion to check two rational endpoint bounds and compute the pairing of their mixtures with a fixed input tensor and output separator. An exact polynomial identity certifies a negative sign on the entire interval [1/3, 5/6]. The program also checks the coordinate retractions in spatial dimensions 2 through 5. The paper supplies the endpoint criterion, and the general dimension extension uses positivity of coordinate restriction and zero padding. [arXiv:2506.14480v1](https://arxiv.org/abs/2506.14480v1)
