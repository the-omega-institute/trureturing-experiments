# Strong q-Laguerre–Pólya class at the zero endpoint

The source is D. K. Dimitrov and B. Shapiro, *Weak and strong q-analogs of
the Laguerre–Pólya class*, arXiv:2606.17864v1, Definition 1.2, Lemma 3.1,
Proposition 3.6 and the zero-endpoint branch of Question 6.5.

Run `python3 check.py` with Python 3.12 or newer. The program uses only
integers and `fractions.Fraction`; it has no external dependencies.
The recorded output is in `output.txt`.

The fidelity checks reproduce the source's multiplication obstruction
at q = 1/2 and verify its reciprocal-tail estimate on the declared finite
grid. The endpoint checks enumerate every multiset of roots in
{-4, -3, -2, -1, 0, 1, 2, 3, 4} through degree six. They verify the
polynomial predictions: at most one nonzero root on each half-axis,
unrestricted multiplicity at zero, exact reconstruction as
c z^m (1+a z)(1-b z), and the corresponding one-sign restriction.
Every excluded multiset has an exact positive rational separation
parameter witnessing its exclusion. Further checks include opposite-sign
roots and origin multiplicities through twenty.

These are finite algebraic checks. They do not establish closedness of
the mesh under analytic limits, the absence of exponential factors, or
the complete Painlevé–Kuratowski limit. Those require the general
analytic argument specified in the preregistration.
