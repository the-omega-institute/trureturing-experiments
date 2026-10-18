# OEIS A394407: every prime other than 3 is a term

A394407 lists the positive k such that no smaller positive m has phi(m)/tau(m) = phi(k)/tau(k). The entry conjectures that every prime p ≠ 3 is a term.

`check.py` (Python standard library) sieves phi and tau below N and checks three things:
- the first 64 published terms;
- every prime p ≠ 3 below N is a term, and 3 is not;
- the inequality phi(m)/tau(m) ≤ (m − 1)/2 for 2 ≤ m < N, which the proof uses.

Run it with `python3 check.py 100000`. Exit 0 and `RESULT PASS` mean all three checks hold. `run.log` is the output for N = 100000.
