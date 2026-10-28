# Chu Conjecture 1.3: the largest MSTD subset of {0,...,n-1} has n - 7 elements

Companion program for H. V. Chu, *The union of two arithmetic progressions with the same common difference is not sum-dominant*, Integers 20 (2020), #A87, Conjecture 1.3 (p. 3); arXiv:1906.03793.

    python3 chu-c13.py    # about 6 minutes; exit 0 iff all checks hold

Checks (exact bit-mask arithmetic):
1. no subset of {0,...,n-1} with at least n - 6 elements is MSTD, exhaustively for 16 <= n <= 40;
2. A_n = {0,...,n-1} minus {3,5,6,n-7,n-6,n-5,n-3} has |A_n+A_n| = 2n-2 > 2n-3 = |A_n-A_n| for 16 <= n <= 40;
3. the six-point lemma |P|+|Q| <= 2|C| for all 1024 pairs X, Y of subsets of {0,...,5} containing 0.
