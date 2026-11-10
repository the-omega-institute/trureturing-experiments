# Blitvić–Steingrímsson Conjecture 5: descent conventions

Source: N. Blitvić and E. Steingrímsson, *Permutations, Moments, Measures*, arXiv:2001.00280v4 (7 October 2020), Sections 3.5–3.6, Definitions 4–8 and Conjecture 5; *Transactions of the AMS* 374 (2021), DOI 10.1090/tran/8330. The arXiv TeX was read. The journal PDF could not be retrieved.

Run `python3 check.py`. It uses only Python's standard library and exact integers. Its captured stdout is `output.txt`; a successful run ends with `ALL_CHECKS_OK`.

The fidelity checks reproduce Proposition 2(1), the counting recurrence for k-arrangements, for 0 ≤ n ≤ 6 and 0 ≤ k ≤ 3, and the explicit 4-arrangement representation example. Both word representations are enumerated independently of colored permutations and checked against their representation maps.

For n = k = 2, the literal color formula with the terminal comparison gives 4t + t². Relabeling fixed-point colors into 0,…,k−1 gives 1 + 4t. The independent permutation and derangement forms both give 3 + 2t. Thus the search candidate's literal mismatch is reproducible.

For ordinary descents comparing only positions i < n, all three tested color labelings (literal, subtract one, and reduction modulo k) give the same polynomial as both word forms for every 0 ≤ n ≤ 6 and 1 ≤ k ≤ 3: 21 parameter pairs, including 18 with n > 0. In particular, deleting the terminal comparison removes the proposed five-object obstruction. The source defines fixed-point colors in [k] but its colored-permutation alphabet is 0,…,k−1, which is also inconsistent without a convention repair.

The candidate is withdrawn under the prescribed rule excluding convention-artifact refutations. The computation establishes finite agreement of the internal-descent variants, not a general identity or the authors' intended convention. No Lean proof, kernel verification or mathematical settlement is claimed.
