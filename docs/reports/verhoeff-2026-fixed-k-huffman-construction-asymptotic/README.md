# Fixed-threshold Huffman construction counts

`check.py` computes the raw symbol counts of Tom Verhoeff's balanced construction
in arXiv:2605.27617v1, Sections 2, 3.1 and 3.4. Section 8, first question,
leaves their fixed-threshold asymptotic open. These counts precede free reduction;
they are not minimum lengths of picture-hanging words.

Run with Python 3 and no external dependencies:

```sh
python3 check.py > output.txt
```

The check uses integers and `fractions.Fraction` throughout. It reproduces the
source's values a[2,2] = 32 and a[2,10] = 157048832, and its formula
a[2,i] = (8/3)6^i - 4·4^i for 2 ≤ i ≤ 160. Greedy merging uses 2(a+b).
It agrees with independent exhaustive optimization over all root bipartitions
on 1,286 multisets with at most eight leaves and weights in {1,2,3,4,8}.
The root-optimization recurrence independently reproduces all 205 counts for
1 ≤ k ≤ 5 and 0 ≤ i ≤ 40, including the auxiliary zero extension outside k ≤ 2^i.

The proposed answer is a[k,i] ~ C_k 6^i i^(k-2), where
C_k = 4^k / (9(k-2)!). The output contains every source-defined count and exact
ratio for k = 3, 4, 5 and i ≤ 40. Displayed decimals are floors of exact rational
numbers at six places, computed by integer arithmetic.

| k | i = 10 | i = 20 | i = 40 | i = 80 | i = 160 |
|---|---|---|---|---|---|
| 3 | 0.592889 | 0.793569 | 0.896759 | 0.948379 | 0.974189 |
| 4 | 0.349454 | 0.631589 | 0.804659 | 0.899544 | 0.949075 |
| 5 | 0.185771 | 0.513442 | 0.731684 | 0.858857 | 0.927591 |

The exact k = 3 identity is checked for 6 ≤ i ≤ 160. The proposed triangular
increment and the resulting form 6^i P_k(i) + b_k 4^i are checked for k ≤ 12
and 40 ≤ i ≤ 160. Polynomial antidifferences are calculated over the rationals;
their integration constants are matched at i = 40 and checked at the other
indices. Their leading coefficients agree with C_k. All assertions pass and
the output ends with `ALL_CHECKS_OK`.

These are finite checks. They establish neither eventual stabilization for all
indices nor the limit for every fixed threshold. Those require the general
weighted-tree argument proposed in the preregistration; no Lean proof is part
of this report.

The program is independently written from the cited recurrence and copies no
third-party implementation. Source reference: https://arxiv.org/abs/2605.27617v1.
