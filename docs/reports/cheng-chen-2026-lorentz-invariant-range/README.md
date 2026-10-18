# Cheng–Chen Conjecture 24: the range of the Lorentz invariant for odd n

Check for the open-problem preregistration on Conjecture 24 of Y. Cheng and L. Chen, arXiv:2605.19530v2: for odd n, the Lorentz invariant I_ρ = Tr(ρ^T ε_n ρ ε_n) of n-qubit states ranges over [−1/2, 0].

Run: `python3 check.py 20000 1` (Python 3 with numpy).

For n = 1, 3, 5 the program checks:
- ε_n^T = −ε_n and ε_n² = −Id;
- 20000 random states of every rank have I real and in [−1/2, 0];
- the endpoints: a pure state gives 0, and (|v⟩⟨v| + |Θv⟩⟨Θv|)/2 with Θv = ε_n v̄ gives −1/2;
- the segment between these two states decreases monotonically from 0 to −1/2.

It also evaluates the endpoints exactly at n = 1, and shows that for even n the range differs (a Bell state at n = 2 gives I = 1).

`run.log` holds the output of the run cited by the preregistration.
