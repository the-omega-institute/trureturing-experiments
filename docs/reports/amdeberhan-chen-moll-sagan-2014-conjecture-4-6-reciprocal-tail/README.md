# Amdeberhan–Chen–Moll–Sagan, Conjecture 4.6: reciprocal tails of {rk}_{s,−t}

Check for the open-problem preregistration proving the first clause of Conjecture 4.6 of T. Amdeberhan, X. Chen, V. H. Moll and B. E. Sagan, "Generalized Fibonacci polynomials and Fibonomial coefficients", Annals of Combinatorics 18 (2014), 541–562 (arXiv:1306.6511v2).

Run: `python3 check.py` (Python 3 with mpmath; 260-digit arithmetic, tails truncated below 10^−240).

With {0} = 0, {1} = 1, {j} = s{j−1} − t{j−2}, the program checks:
- for 2 ≤ s ≤ 12, 1 ≤ t < s, (s, t) ≠ (2, 1), r ∈ {1, 2, 3, 5} and n ∈ {1, 2, 3, 4, 7} (1300 cases): ⌊(Σ_{k≥n} 1/{rk})^{−1}⌋ = {rn} − {r(n−1)} − 1, and the proof's bounds 1/b_n < S_n ≤ 1/(b_n − 1) with b_n = {rn} − {r(n−1)};
- exactly, for 2 ≤ s ≤ 9 and r ≤ 4, the Cassini identity a_k² − a_{k−1}a_{k+1} = c²Q^{k−1} and the bound a_{k+1} − a_{k−1} ≥ c²Q^{k−1} + 1 for a_k = {rk}, c = {r}, Q = t^r;
- the excluded pair (s, t) = (2, 1) gives {n} = n.

The tail precision matters: with 80 digits and a 10^−70 cutoff one case (s, t, r, n) = (12, 7, 5, 7) is misjudged, because 1/S_n lies within a relative 10^−37 of the integer b_n.

`run.log` holds the output of the run cited by the preregistration.
