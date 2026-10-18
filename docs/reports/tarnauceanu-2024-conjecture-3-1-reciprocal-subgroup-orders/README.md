# Tărnăuceanu, Conjecture 3.1: the infimum of f over finite non-cyclic groups is 1

Check for the open-problem preregistration proving Conjecture 3.1 of M. Tărnăuceanu, "On the sum of inverses of subgroup orders in finite groups", Scientific Studies and Research, Series Mathematics and Informatics 34 (2024), no. 1, 115–120.

Run: `python3 check.py` (Python 3, standard library; exact rationals).

For f(G) = Σ_{H ≤ G} 1/|H| and the split groups G = C_p ⋊ C_{q^s} (q | p − 1, the generator acting by an automorphism of order q), the program checks:
- by brute-force subgroup enumeration for five small parameter triples, that G has 2s + p + 1 subgroups, is nonabelian, and that f(G) equals Σ_{j=0}^{s} 1/(p q^j) + Σ_{j=0}^{s−1} 1/q^j + p/q^s;
- for q ∈ {3, 5, 11, 23, 47, 101}, with p the least prime ≡ 1 (mod q) and s the least exponent with p/q^s < 1/q, that 1 < f(G) < 1 + 3/(q − 1).

`run.log` holds the output of the run cited by the preregistration.
