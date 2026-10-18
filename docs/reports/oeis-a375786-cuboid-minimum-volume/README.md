# OEIS A375786: the cuboid with the smallest edge has the minimum volume

A375786(n) is the minimum volume of an integer cuboid x ≤ y ≤ z whose surface equals that of the cube of edge n, that is, xy + yz + zx = 3n². The entry conjectures that the cuboid with the smallest edge, with ties broken by the smallest second edge, always has the minimum volume.

`check.py` (Python standard library) enumerates every solution for each n ≤ N through the divisors of x² + 3n². It checks four things:
- the conjecture itself;
- the first 36 published terms;
- the inequality 2xyz ≥ x(S − x²) on every solution, with S = 3n²;
- that the least edge is 1, or else 2, in which case 1 + 3n² is prime, n is even and the second edge is 2.

Run it with `python3 check.py 1500`. Exit 0 and `RESULT PASS` mean every check holds. `run.log` is the output for N = 1500 (19 s).
