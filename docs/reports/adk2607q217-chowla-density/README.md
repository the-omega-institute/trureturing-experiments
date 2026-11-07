# Chowla extremal-equality density: exact checks

Source: Aliabadi, Driscoll, Krop, Sirkovic, Sullivan and Vedadi, *Extremal Chowla Sets and Their Linear Analogues*, arXiv:2607.24847v2, Definition 2.1, Theorem 2.4, Proposition 2.6 and Question 2.17.

Run `python3 check.py`. The script uses Python's standard library, integers and `fractions.Fraction`; it uses no floating point. Successful output is stored in `output.txt` and ends with `ALL_CHECKS_OK`.

For every modulus from 2 through 4,000, orders of all residues determine the maximum Chowla-set size. This is checked against the source's divisor formula and its equality criterion. Literal subset enumeration supplies an independent oracle for moduli 2 through 12.

For prime cutoffs 7, 11 and 13, a complete residue-class count modulo the primorial equals the finite Bernoulli prime-pattern mixture. Finite-pattern products are also checked for distinctness. Exact empirical ratios are given at 4,000, 20,000 and 200,000. The proposed rational enclosures for the limiting density at parameters 2, 3 and 4 contain the last finite ratio. This containment supplies a consistency check; a finite ratio has no certified error estimate to the limit here.

The product of the primes from 3 through 23 supplies an exact strict-inequality witness for a follow-up upper bound. Finite computations establish neither existence of the natural density nor atomlessness of an infinite distribution; those require the general argument in the preregistration.
