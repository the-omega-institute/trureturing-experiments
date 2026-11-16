# `liu2503s6` exact check

`check.py` enumerates all three-voter Luce ranking triples for `n = 4` with exact `Fraction` arithmetic and weights `(1/6, 1, 1, 1)`. It checks the Luce probability normalization and the exponential-race fidelity value `P(0 beats 1) = 1/7`, computes the exact probability that alternative `0` is a Condorcet loser and the resulting Smith-set sizes, and verifies the fourth-power form of the proposed lower bound `(1/2)(6n)^(-1/4)` without floating point. The finite computation is a sanity check; the preregistered settlement uses the all-`n` analytic lower bound.
