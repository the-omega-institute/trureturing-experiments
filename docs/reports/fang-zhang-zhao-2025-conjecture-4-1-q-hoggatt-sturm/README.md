# Fang–Zhang–Zhao, Conjecture 4.1: the q,t-Hoggatt polynomials form a Sturm sequence

Check for the open-problem preregistration proving Conjecture 4.1 of H. Fang, C. X. T. Zhang and J. J. Y. Zhao, "Analytic properties arising from the Baxter numbers", Journal of Mathematical Analysis and Applications 561 (2026), 130615 (arXiv:2505.05873v1).

Run: `python3 check.py` (Python 3 with sympy; exact rationals and exact real-root isolation; about one minute).

With M(k, l, m; q) = q^{m·binom(k+1, 2)} Σ_{π ∈ PP_{k,l,m}} q^{|π|} and H_n^{[m]}(q, t) = Σ_{k+l=n−1} M(k, l, m; q) t^k, the program checks:
- MacMahon's box product against brute-force enumeration of plane partitions (q = 2, 3; m = 2, 3; k, l ≤ 3);
- Lemma 1 of the proof: H_n^{[m]} = F_{n+m−2,m} / ∏_{j<m} C_q(n+m−2, j) (q = 2, 3; m = 2, 3; n ≤ 7);
- Lemma 5: D_q(F_{N,1} ∗ F_{N−1,m−1})(t) = [N]_q q^m F_{N−1,m}(q^m t) (q = 2, 3; m = 2, 3, 4; N ≤ 9);
- strict interlacing of consecutive H_n^{[m]}, n ≤ 10, all zeros negative, for q ∈ {3/2, 2, 3} and m ∈ {2, 3, 4}.

`run.log` holds the output of the run cited by the preregistration.
