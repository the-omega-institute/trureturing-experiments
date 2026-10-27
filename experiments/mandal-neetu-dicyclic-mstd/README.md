# Mandal–Neetu Conjecture 6.2: many MSTD sets in dicyclic groups

Companion program for S. Mandal and Neetu, *Exact and Asymptotic Counts of MSTD, MDTS, and Balanced Sets in Dicyclic Groups*, arXiv:2602.09073v1, §6, Conjecture 6.2.

    python3 mn-c62.py    # checks |A_t A_t| = 4r, |A_t A_t^{-1}| = 4r-1 for A_t = {a,...,a^r} u {a^t b}, r = m-1, m = 2..7, n = 2m-1..2m+6, and brute-forces S_n(m) >= 2n for five small (m, n); exit 0 iff all hold

Exact integer computation in the model (i,e) <-> a^i b^e of Dic_{4n} = <a,b | a^{2n}=1, b^2=a^n, b^{-1}ab=a^{-1}>.
