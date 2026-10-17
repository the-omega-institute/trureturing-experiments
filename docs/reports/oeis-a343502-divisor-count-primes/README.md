# OEIS A343502: tau(tau(k)) and tau(k+1) both prime

Check for the open-problem preregistration on Dequatre's conjecture in OEIS A343502: all members except 2 and 3 are composite.

Run: `python3 check.py 2000000` (Python 3, standard library only).

The program sieves the number of divisors up to N + 1 and lists the members k ≤ N. It checks:
- the published data;
- that the only prime members are 2 and 3;
- that every odd member k has k + 1 = 2^(t−1) with t prime, the shape used in the proof.

`run.log` holds the output of the run cited by the preregistration.
