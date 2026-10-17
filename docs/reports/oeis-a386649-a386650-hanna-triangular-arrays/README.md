# OEIS A386649 and A386650: products of trinomial coefficients and Hardin's triangular arrays

A386649(n) is the product of the central trinomial coefficients [x^k](1 + x + x²)^k for k ≤ n. A386650(n) is the same product for the quadrinomial coefficients [x^k](1 + x + x² + x³)^k. Paul D. Hanna conjectures:
- A386649(n) = A214589(n) − 2;
- 2·A386650(n) = A214590(n) − 2.

Here A214589 and A214590 are R. H. Hardin's counts of n × n × n triangular arrays with entries 0..2 (resp. 0..3) whose rows, of lengths 1..n, all have the same average.

`check.py` (Python standard library) does three things:
- counts the arrays by brute force for small n;
- counts them for every published n by summing, over the averages 0..d that the one-entry first row allows, the product of the row counts;
- compares both counts with the data and with the product formulas.

Run it with `python3 check.py`. Exit 0 and `RESULT PASS` mean everything agrees. `run.log` is the output.
