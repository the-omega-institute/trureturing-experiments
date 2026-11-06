# arXiv:2610.18636v2, Problem 21

This report is a preregistration-stage exact arithmetic check for Jovanović, Sedlar and Škrekovski, *Disproving two conjectures on the cospectrality of complete tripartite graphs*, §5, Problem 21.

`check.py` uses Python integers, `fractions.Fraction`, and SymPy exact algebra only. It checks the cubic endpoint identities used for the Hₙ spectrum, derives the published Hₙ/Kₙ,ₙ,ₙ spectral-distance identities symbolically, verifies the n=3 Hₙ characteristic polynomial and rank, checks the connected complete-graph rank-six-or-more family, and verifies the two nominated integer certificate numerators 750 and 1366 over denominator 1024. The saved output is `check.out` and ends with `ALL_CHECKS_OK`.

These are fidelity and certificate-arithmetic checks. They do not constitute the general rank/interlacing proof, which is the post-registration probe.
