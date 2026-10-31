# Au G₂ Witten-zeta nonvanishing: exact checks

Sources: Kam Cheong Au, *On single-variable Witten zeta functions of rank two and three*, [arXiv:2412.17196v4](https://arxiv.org/abs/2412.17196v4), Example 5.5 and Conjecture 5.6; [v1](https://arxiv.org/abs/2412.17196v1), section 4, the displayed recurrence for g_n.

Run `python3 check.py`. The program uses only integers and `fractions.Fraction`, needs no third-party packages, and exits nonzero on a failed assertion. The stored result is `output.txt`; its final line is `ALL_CHECKS_OK`.

The fidelity checks reproduce the quartic reversal identity, the vanishing coefficient classes in Example 5.5, and 41 instances of the original published g_n recurrence, using generalized-binomial convolution independently of that recurrence. The settlement checks compare the literal terminating Gauss sum with the proposed cleared-denominator integer N and its odd parity on 1,000 pairs (0 ≤ m ≤ 24; odd -49 ≤ k ≤ 49, 5 ∤ k), and on every admissible conjecture parameter -299 ≤ k ≤ 1 (121 values).

These finite exact checks are diagnostics for a proposed uniform proof, not a proof of the universal conjecture. They include m = 0. No Lean proof or analytic continuation calculation is included.
