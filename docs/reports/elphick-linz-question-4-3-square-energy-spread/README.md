# Elphick–Linz Question 4.3: negative quadratic squared spread

Source: Clive Elphick and William Linz, *Symmetry and asymmetry between positive and negative square energies of graphs*, Electronic Journal of Linear Algebra 40 (2024), 418–432, [DOI:10.13001/ela.2024.8447](https://doi.org/10.13001/ela.2024.8447); [arXiv:2311.11530v3](https://arxiv.org/abs/2311.11530v3), Question 2. The journal numbers it Question 4.3.

Run `python3 docs/reports/elphick-linz-question-4-3-square-energy-spread/check.py` from the repository root. The program uses Python 3 and SymPy; the stored output uses Python 3.14.8 and SymPy 1.14.0. All evidence uses integers and exact rational arithmetic.

The fidelity check independently computes the complete-graph characteristic polynomials for orders 1 through 12 and reproduces the source's squared-spread identity `(n−1)(n−2)`. The candidate check constructs the disjointness graph on the ten two-subsets of a five-element set, verifies connectivity and `A²=2I−A+J`, and computes its complete characteristic polynomial. For every equal independent blowup of size 1 through 8, it directly constructs the adjacency matrix, computes all characteristic roots with multiplicity, and verifies connectivity and square energies `s⁺=14t²`, `s⁻=16t²`. Thus every checked squared spread is `−2t²=−n²/50`.

These finite checks support the proposed construction. They do not prove the assertion for unbounded orders or establish literature novelty. The universal argument uses the full decomposition into vectors constant on each fibre and vectors with zero sum on each fibre. `output.txt` stores the successful run and ends in `ALL_CHECKS_OK`.
