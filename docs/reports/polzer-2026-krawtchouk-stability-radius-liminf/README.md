# Krawtchouk stability-radius diagnostics

The target is the positive-liminf subquestion in Leslie P. Polzer, “A log-concavity criterion for finite nonnegative linearization”, arXiv:2611.01561v1, Section 6. The optimal exponential rate is outside this target. The reference margin is Lemma 3.1 of the companion manuscript “Exact and perturbative criteria for finite nonnegative linearization”, release v1.0.7 of skypher/szkraw.

Run from the repository root:

```sh
python3 docs/reports/polzer-2026-krawtchouk-stability-radius-liminf/check.py
```

Python 3 standard library only. All evidence uses integers and fractions.Fraction. The stored output is output.txt; the expected last line is ALL_CHECKS_OK.

For each N=1,...,8 at r=1/2 and r=1/3, an independent subset-product formula agrees with every matrix linearization coefficient. This reproduces the published corner identity g(N,N;N)=(1-r)^N, row sums, triangle zeros and the reference margin (min(r,1-r)/4)^N.

The program checks the proposed radius eta=(1/6)(min(r,1-r)/72)^N on 480 sample arrays: the reference, every one-row hexagon vertex, both terminal-row directions, and six coupled choices. It checks the exact subtracted matrix recurrence, its scaled norm inequality, the factorial lower bound for the leading-coefficient product, every coefficient error bound and the positive half-margin. Repeated arrays at N=1 are included in the sample count. These are finite diagnostics, not an exhaustive test of a real coefficient ball and not a proof of the unbounded assertion. No Lean proof is included.

Producer: independent codex-cli research seat, flight 85a1c20d458e04c23a439a8e, repo-prior-exposed. The program is an independently written implementation of the displayed formulas, with no copied source code. The paper and companion release supply the cited mathematical identities and margins.

check.py SHA-256: 472ea4cb71514594b22430cc34019388db5a82107887425071f3fb200caa4962
