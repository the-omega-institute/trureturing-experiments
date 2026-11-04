# Gaussian quotient monotonicity: Problem 7.2

Exact arithmetic certificates for the final yes/no question in Madjid Mirzavaziri and Daniel Yaqubi, *Graph binomial coefficients and a positivity conjecture of Guo and Krattenthaler*, arXiv:2610.19249v1, Section 7, Problem 7.2. The source is https://arxiv.org/abs/2610.19249v1.

Run `python3 check.py`. Python 3's standard library is sufficient; all evidence uses integers. `check-output.txt` is the resulting output and ends with `ALL_CHECKS_OK`.

Fidelity checks reproduce Theorem 4.3(iii): the negative positions of the source's P^(4) are 1, 5, 59, 63, and those of P^(6) are 1, 7, 147, 153, all with coefficient -1. The explicit source identity P^(1)=1+q²+q⁴ is also checked. Gaussian products agree with an independent q-Pascal implementation in the fidelity cases and the small counterexample.

The counterexample parameters are (N,k,a,b)=(30+420m,3,7,10). For m=0, the program independently enumerates all 4,060 three-element subsets of {0,...,29}, weighted by sum(S)-3. Exact polynomial multiplication checks U[7]=G and P[10]=U. U increases weakly from degree 1 to its middle, while P has negative coefficients at degrees 1, 11, 55, 65.

The program verifies the two 42-entry residue tables, the degree-32 polynomial E=[42]/([2][3][7]) with E(1)=1, and the difference recurrence through 50 periods. It also builds the full Gaussian quotient coefficients for m=0,...,20 and m=100, checking divisibility, degrees, reciprocity, first-half monotonicity, the difference formula and p_11=-1 in every case. These finite checks support the separate ordinary argument for all m; they do not establish unbounded quantification or claim a Lean proof. The broader characterization in Problem 7.2 and Conjecture 4.4 are outside the result.
