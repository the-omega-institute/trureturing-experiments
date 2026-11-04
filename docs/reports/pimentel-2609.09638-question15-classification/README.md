# Mean-distance minimizers on the sphere

Exact helper checks for J. V. M. Pimentel, *The mean distance to a simple closed curve on the sphere*, arXiv:2609.09638v1, Question 15 (classification clause).

Run `python3 check.py` with Python 3 and SymPy. The recorded run uses SymPy 1.14.0; its output is in [output.txt](output.txt).

The program reproduces the source's six-arc construction: sphere membership, unit speed, matching points and tangents, and total length (4\pi). It checks every endpoint reconnection graph for (2\le n\le32), the latitude-length identity in exact cyclotomic arithmetic for (2\le n\le12), the source's mean-distance values for (n=2,3,4) using its published tube formula, and a rational instance of the proposed rigidity lemma's perturbation construction. Assertions use exact arithmetic. Decimal values for (n=2,3) display exact expressions and supply no additional evidence.

These are fidelity and helper checks. They do not verify the global tube formula, the area formula for rectifiable curves, the published thick-rope classification, or the proposed universal classification of minimizers. No Lean proof or probe is included.
