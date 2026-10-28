# Finite preparation certification

`check.py` uses Python 3 standard-library rational arithmetic. Run it with
`python3 check.py`; `check.out` contains its output. It exits successfully only
after printing `ALL_CHECKS_OK`.

Source: Yujie Zhang, Yìlè Yīng and David Schmid, *Quantifiers and witnesses for
the nonclassicality of measurements and of states*, arXiv:2504.02944v3,
Quantum 10, 2180, DOI 10.22331/q-2026-07-30-2180.

The fidelity check reproduces the entire rational probability table and the
noncontextual ontological model in Appendix E. It checks bases of all linear
preparation and measurement identities, not only selected identities. The
diagonal witness value is zero, even though the displayed finite experiment
has a model. Thus this check retains the source's distinction between a
theory-dependent SDP witness and a noncontextuality inequality.

The independent checks use the noisy four-outcome qubit measurement
`(I + eta X)/4`, `(I - eta X)/4`, `(I + eta Z)/4`, `(I - eta Z)/4`
and four states `(I + u a X + v b Z)/2`, where `u,v` are signs. The inequality
`W = sum(u*(p_x+ - p_x-) + v*(p_z+ - p_z-)) <= 2`
is checked on all four normalized preparation-cone rays paired with all four
response vertices. Its quantum value is exactly `2*eta*(a+b)`.

All 1,677 cases with `a,b,eta` reduced fractions in `[0,1]` of denominator at
most 6 and `a*a+b*b <= 1` pass: 1,480 nonzero-noise cases have explicit models,
68 violate the inequality, and 129 have zero noise. The 21 equality cases
are included in the 1,480 models. Extra identities at degenerate parameters
are checked by exact row reduction. Separately, each of the 95 distinct
rational noise parameters of denominator at most 32 satisfying `2*eta*eta > 1`
has a witness consisting of four rational pure states. A final check compresses
25 responses to four common vertices and verifies preservation of statistics,
normalization, positivity and the preparation relation.

These are finite checks of the proposed mechanism. They establish neither
Conjecture 1 for all measurements nor the proposed compactness argument or
preparation bound. The planned positive settlement is a general argument;
these computations are external supporting experiments.
