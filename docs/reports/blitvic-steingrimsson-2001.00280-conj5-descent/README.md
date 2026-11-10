# Descent conventions for Blitvić–Steingrímsson Conjecture 5

Source: *Permutations, Moments, Measures*, arXiv:2001.00280v4, Sections 3.5–3.6, Definitions 4–8 and Conjecture 5; Trans. AMS 374 (2021), DOI 10.1090/tran/8330. The journal text was not read.

The published definitions have an unresolved convention conflict. Fixed-point colors in the encoding range over 1,…,k, while the colored-permutation alphabet is 0,…,k−1. Definition 5 explicitly counts the terminal comparison with (n+1,0). Omitting that comparison gives a different statistic. This report supports **abstention from filing either a proof or a refutation of the published clause** until the source convention is resolved.

Run with Python 3.12 or newer:

```sh
python3 check.py
```

The script uses only the Python standard library and integer arithmetic. Its complete output is in `output.txt`, ending in `ALL_CHECKS_OK`.

The fidelity checks reproduce the source's arrangement-count recurrence and permanent formula for 0≤n≤9, 0≤k≤4, and its explicit 4-arrangement representation example. They also reproduce the type-B size-one polynomial 1+t with the terminal convention; internal comparisons alone give 2.

For 0≤n≤9 and 1≤k≤4 the program computes every coefficient of both word-form descent distributions and of the colored distributions with and without the terminal comparison. It tests the literal positive colors and the two zero-based relabelings φ−1 and φ mod k. The latter have the same full color alphabet. All three internal distributions equal both word distributions. For n=k=2 the word polynomial is 3+2t; literal terminal colors give 4t+t², and either zero-based terminal labeling gives 1+4t.

These are exact finite counts. Colored counting enumerates every permutation and groups it by its fixed-point mask and ordinary descent mask, then sums every coloring with a transfer recurrence. The word computation is independent: it tracks unused positive letters, remaining negative slots and the previous letter, excluding fixed points of the positive subword for derangement forms. Explicit enumeration through n=5, k=4 checks the transfer computations, the independently generated word sets, the source representation maps and the proposed internal-statistic bijections. No sampling or numerical approximation is used.

A general bijection for the **modified internal statistic** is also checked object by object through n=5. Reverse-complement the underlying permutation, mirror the fixed-point colors, replace positively colored fixed points by the negative of their colors, and standardize the remaining positive letters increasingly. With literal positive colors this lands in derangement forms; with zero-based colors it lands in permutation forms. Each adjacent comparison is preserved, as a comparison at the reversed position. This does not preserve the extra terminal term and is not a proof of the source clause. No Lean proof was written or compiled.
