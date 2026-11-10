# Exact checks for Conjecture 2.4

Source: Siva Athreya, Satyaki Mukherjee and Soumendu Sundar Mukherjee, *On a generalisation of the coupon collector problem*, arXiv:2304.01145v2, Conjecture 2.4.

`check.py` uses only integers, `fractions.Fraction`, and exact SymPy rational algebra. It independently solves the finite birth-death first-step equations for a fixed target s-subset and verifies the recurrence and its symbolic leading coefficient for the requested `(r,s)` pairs through `n = 12`. For small `(n,r,s)` it computes the complete Markov cover expectation by exact rational layer-by-layer linear systems and checks the tied-discovery record-cover squeeze. It also checks the source's Johnson-graph degree statement for the proved `r=s` case. The output is in `check.out`; the last line is `ALL_CHECKS_OK`.

The exact checks are consistency and fidelity checks only. They do not constitute a proof of the asymptotic conjecture.
