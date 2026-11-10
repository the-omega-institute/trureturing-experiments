# Rational normality-degree realization checks

Source: Marius Tărnăuceanu, *Normality degrees of finite groups*, arXiv:1312.1427v1, §1 and Conjecture 4.4; Carpathian Journal of Mathematics 33 (2017), 115–126, Conjecture 4.1. DOI: 10.37193/CJM.2017.01.12.

Run `python3 check.py`. Python 3.12 or newer and the standard library suffice. All evidence uses integers, sets and `fractions.Fraction`; no floating point is used. The stored result is `output.txt`.

The fidelity check independently enumerates the six subgroups of S3, finds three normal subgroups and reproduces the source's value 1/2. The general checker constructs the multiplication table of F_p semidirect C_(2^B), exhausts its entire subgroup lattice by adjoining generators and tests normality using every ambient conjugator. This enumeration does not use the proposed classification. It then compares the output with |N| = 2B − t + 2 and |L| = 2B − t + 2 + pt, and checks the predicted classification separately at each quotient-subgroup level.

The direct cases use p in {3,5,7,11,13,17}, every t in {1,2,3,4} for which the field has an element of order 2^t, and every t ≤ B ≤ 4, choosing the least suitable action element. One additional direct case is (p,B,t,r) = (5,5,2,2), the full realization recipe for 1/2. Larger parameters are checked only by exact arithmetic: for every 1 ≤ u < v ≤ 10, find a prime p congruent to 1 modulo 2^(2(v−u)), construct an element of that exact order, and check B = pu + v − u − 1 and the predicted count identities. These arithmetic checks do not enumerate the corresponding large groups.

The checks support a proposed general proof, and do not establish the universal conjecture. No Lean proof is supplied. The programs are independently authored from the mathematical definitions; no third-party program source is copied.
