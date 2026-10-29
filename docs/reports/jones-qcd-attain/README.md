# Jones quantum-channel fidelity attainment

`check.py` checks the candidate negative answer to the positive-distance attainment question in Morgan Jones, *Continuity of certain quantum channel divergences*, arXiv:2610.16270v2, Section 13(1), concerning Theorem 4.6 only.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The stored `output.txt` ends in `ALL_CHECKS_OK`. All calculations use exact integers, rationals and algebraic square roots, with no floating point or tolerance.

The source fidelity check reproduces Proposition 3.2 for orthogonal pure-state replacement channels: the unequal splitting has fidelity zero and the equal splitting fidelity one. The candidate checks compute the unnormalised Choi matrices from matrix units, verify Kraus representations and trace preservation for four nonzero channels at five rational parameters, and verify the refined quadratic and square-root estimates on 196 ordered effect pairs, including 50 noncommuting pairs. A free noncommutative algebra calculation checks the anticommutator identity.

These are finite exact checks and symbolic identities. They do not establish the universal nonattainment theorem, the PSD geometric-mean limit construction, or its Lean formalization. The proposed proof requires those separate arguments. Source reference: https://arxiv.org/abs/2610.16270v2.
