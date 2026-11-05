# Alpha-triangle integrality: exact finite checks

`check.py` checks Suresh Govindarajan, *Notes on higher-dimensional partitions*, arXiv:1203.4419v1, Section 2.5 and Conjecture 3.1 (Section 3.1). Run with Python 3 and no third-party packages:

```sh
python3 check.py
```

The program independently enumerates all labelled finite coordinatewise lower ideals obtained by adding at most four cells to μ_r = {0,e₁,…,e_r}, for every 0 ≤ r ≤ 8. Immediate-predecessor tests decide valid additions. Filtering by full support reproduces the source's entire C-matrix rows m = 0,…,4. These source rows are fidelity assertions, not enumeration inputs.

For every strict skew diagram in those rows, coordinate permutations independently check degree-two coverage, the action kernel's injection into a product of two-element permutation groups, stabilizer-order divisibility by 2^m m!, and orbitwise denominator cancellation. Exact `Fraction` polynomial arithmetic verifies monomial-coefficient integrality and the binomial transform against the independently enumerated counts for every 0 ≤ m ≤ 4 and 0 ≤ r ≤ 8. The sharp disjoint-edge count is also checked in every row.

`output.txt` is the program's complete output. Its last line is `ALL_CHECKS_OK`. The finite checks support the proposed proof route; they do not prove the unbounded conjecture. No Lean proof or probe is included.
