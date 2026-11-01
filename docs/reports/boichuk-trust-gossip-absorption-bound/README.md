# Boichuk trust-gossip absorption bound

Source: Nicholas Boichuk, “Termination of Binary Trust-Gossip Dynamics: A Constructive No-Limit-Cycles Theorem”, arXiv:2607.04023v1, submitted 4 July 2026, Section 10, Question 1. Source: https://arxiv.org/abs/2607.04023v1

Run with Python 3 from this directory:

    python3 check.py

`output.txt` contains the complete deterministic output. All state encodings, reachability calculations, move checks and breadth-first distances use exact integers; random states use seed 260704023 and integer sampling.

The fidelity check constructs the complete state graph for n=3,4 and uses reverse breadth-first search along flat firing moves to find the shortest preparation before absorption or an available progress move. It reproduces the source values 1 and 2. Reachability counts exclude diagonal pairs, and every firing uses three distinct agents.

The construction check exhausts all 4,165 states for n=1 through 4 and adds 780 structured and seeded states at n=5,6,10,20,50. It completes a sink component, aligns its entire ancestor basin, and recurses on the remaining agents. A separate dictionary-based replay checks every move against the original update rule, final absorption, the n(n-1) bound, and the stronger condition that no off-diagonal coordinate changes twice.

These are finite checks supporting a proposed general proof. They do not constitute a Lean proof or establish the universal bound. They give no upper bound on every fair schedule or on expected random absorption time, and make no optimality claim for the bound.
