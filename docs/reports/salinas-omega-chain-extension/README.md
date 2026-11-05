# Omega-chain extension checks

Source: Remi Salinas Schmeis, *General Criterion for a Poset to Admit a
Partition into omega-Chains*, arXiv:2610.29056v1, Question 2.32 and Theorem 2.5.

Run with Python 3, without dependencies:

```sh
python3 docs/reports/salinas-omega-chain-extension/check.py
```

The program uses integer indices and Boolean order relations. It enumerates
every partial order compatible with a fixed linear-extension labeling on
one through five elements: 407 orders, containing every isomorphism type in
that size range.

For source fidelity, it constructs the countable order P + omega and checks
the partition with parts p_i, t_i, t_(n+i), t_(2n+i), and so on. These are
instances of Theorem 2.5. It checks induced embeddings, order axioms, the
displayed chain prefixes and covering/disjointness, and strict successor
formulas beyond the finite window. It does not enumerate all countable posets.

For the proposed negative answer, it checks separated-ray and merged-tail
extensions of connected stars with two through six maxima and ray heights
one through eight. In the separated order, incomparable ray points have no
common upper bound. In the merged order, a strict successor of one member
of an incomparable pair is above the other member. These are the local
relations used by the proposed embedding obstruction.

The additional finite-classification checks cover all 868 choices assigning
each base point to a maximal element above it, test the resulting ray
partitions, and find no induced copy of a finite poset without a greatest
element inside any one maximal lower cone plus a ray (912 components).
For these finite copy tests, a ray can be compressed to n points: an n-point
copy selects at most n ray points and only their relative order matters.
The program also checks 337295 explicit merged-tail escape witnesses and
the proposed greatest-element embedding on the tested prefixes.

These checks are finite exact evidence. They do not establish the infinite
partition or least-extension assertion. The unbounded mathematical argument
and its Lean verification remain separate obligations. The stored output is
in output.txt; its final line is ALL_CHECKS_OK.
