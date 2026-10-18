# Giacomelli, Section 7 (1): is the cone of infinitely log-concave initial data convex?

Check for the open-problem preregistration answering the convexity question of P. Giacomelli, arXiv:2604.14391v1, Section 7, item (1).

Run: `python3 check.py 40` (Python 3, standard library only; exact rational arithmetic).

For the order-3 recurrence a_{n+1} = (n+3)a_n − (3n+2)a_{n−1} + 2n·a_{n−2} (n ≥ 2), the program checks:
- the initial data (1, 1, 1) and (1, 2, 4) generate u_n = 1 and v_n = 2^n;
- every iterate of the log-concave operator on u and v is ≥ 0, since L(u) = L(v) = 0;
- the midpoint (1, 3/2, 5/2) generates (1 + 2^n)/2, whose first iterate is −1/4 at n = 1 and −1/2 at n = 2.

`run.log` holds the output of the run cited by the preregistration.
