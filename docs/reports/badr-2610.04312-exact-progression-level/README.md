# Exact checks for Badr's progression-level question

`check.py` uses Python 3.12+ and SymPy 1.14, with integer and rational arithmetic only. Run `python3 check.py`; `output.txt` contains the output.

The fidelity check reproduces the worked quintic in Section 9 of Eslam Badr, *Cubic Progressions on Smooth Plane Curves*, arXiv:2610.04312v1: its rational point, genus 6, projection degree 4 and four simple zeros of X/Y. It also checks smoothness by exact projective-chart Groebner bases.

The candidate checks verify D/e <= d/2 < d-1 for 4<=d<=100 and 2<=e<=100, and seek an irreducible full-degree fiber in the first 16 terms of each sampled progression for two smooth curve families of degrees 4 through 8. Parameters are a,b in {-2,-1,1,2}, with torsion geometric ratios excluded. A full-degree irreducible fiber obstructs every lower level in each sampled configuration. Both projection centers on and outside the curve are tested.

These finite checks do not establish the universal threshold, Hilbert irreducibility, geometric gonality, or pullback irreducibility. The proposed positive settlement needs a general mathematical argument. The source is https://arxiv.org/abs/2610.04312v1; the verification program is independently authored and copies no third-party code.
