# Exact checks for Clow–van Bommel, Conjecture 8.2

Run `python3 check.py` and compare stdout with `output.txt`. The script uses Python's standard library and integer arithmetic.

For fidelity, it recomputes the eternal parameter in three instances of Theorem 7.2 by enumerating guard multisets, matching guard occurrences across moves, and repeatedly removing configurations that fail to respond to some attack. It then computes static domination by exact set cover on distance balls and checks every local response of the block strategy. The latter checks cover odd and even lengths on three cores and six members of the subcubic counterexample family.
