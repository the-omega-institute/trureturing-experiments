# Exact chirp-constant checks

Source: JiYe Liu, Wenkai Wang, Qiang Tian and Wenjun Wang,
[Cluster-Graph Edit Distance: Optimal Explicit Embeddings, Metric Proxies, and Complexity](https://arxiv.org/abs/2608.17990v2), Proposition 7.1 and Conjecture 7.2.

Run with Python 3 and SymPy 1.14.0:

```sh
python3 docs/reports/liu-wang-cluster-edit-chirp-constant/check.py
```

The stored output is `check-output.txt`. The check uses integers, rational arithmetic and exact radicals throughout. Dyadic levels run from zero through the exact integer ceiling logarithm of n, including all zero padding to the resulting power of two. A sparse sum omits only zero coordinates, so large padding does not require a large array.

The fidelity checks reproduce the single-chirp partition sizes, profile, norm, total variation and every level energy for m = 2, 6, ..., 82. In particular the n = 30 example has level energies (340, 20, 4, 4, 4, 4) and critical energy 343 + (23/2)sqrt(2). Exact symbolic upper and lower bounds recover the source's single-chirp limit 2/3.

The proposed refutation uses the multiset union of the chirps at m = 3k and m = k, with k congruent to 2 modulo 4. Direct partition evaluation checks the proposed size, l1 norm, level energies and tail estimate. Exact symbolic bounds have common limit sqrt(2/5), strictly below 2/3. A separate rational estimate bounds the entire family by 13/20 for k >= 50; its rational upper bound is 1250678/1953125. This estimate supplies an unbounded subsequence of admissible sizes, and the infimum over all distinct partition pairs at each such size is bounded above by the displayed pair's ratio.

Additional checks cover all pairs and triples of chirps from m in {2, 6, ..., 30}, and reciprocal-odd nested families with one through six layers, including common singleton padding across a change of dyadic depth. These additional checks are supplementary; the planned settlement only needs the two-layer family. Neither sharpness of the improved bound nor existence of the true optimal limit is asserted.

Finite checks and exact algebra support a preregistration outline. This directory contains no Lean proof and does not establish a kernel-verified settlement.
