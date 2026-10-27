# Goh Conjecture 15: entropic Glibichuk–Konyagin lemma

Companion program for M. K. Goh, *On an entropic analogue of additive energy*, Essent. Number Theory 5 (2026), no. 2, 243–269, doi:10.2140/ent.2026.5.243, arXiv:2406.18798v4, Conjecture 15.

    python3 goh-c15.py    # exact law of Z = sum_{i<=k} X_i X'_i for mu = (1-1/k) delta_0 + (1/k) U_p; exit 0 iff (89) fails for every tested (k, p)

Tested pairs: (k, p) = (3, 67), (3, 257), (4, 257), (5, 1031). The laws are exact rationals; the entropies (bits) and margins are double precision and are computed observations.
