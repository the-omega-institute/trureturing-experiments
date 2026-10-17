# OEIS A376114: refactorable numbers that are twice a square

A376114 lists the refactorable numbers (τ(k) | k) of the form 2m². Its entry conjectures a(n) = A181795(n)/2. A181795 lists the k whose number of odd divisors is an odd divisor of k and whose number of even divisors is an even divisor of k.

`check.py` (Python standard library) builds both sequences below N from their definitions and checks:
- the conjecture, term by term;
- that every term of A181795 is an even square;
- the proof's criterion. Write m = 2^a·u with u odd, and let c be the odd part of 2a + 2. Then 2m² is refactorable iff 4m² ∈ A181795 iff c·τ(u²) divides u².

Run it with `python3 check.py 2000000`. Exit 0 and `RESULT PASS` mean every check holds. `run.log` is the output (164 terms, 7 s).
