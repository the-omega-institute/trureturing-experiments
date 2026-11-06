# Binary palindromic periodicity enumeration

Source: Gabriele Fici, Jeffrey Shallit, and Jamie Simpson, *Some Remarks on
Palindromic Periodicities*, arXiv:2407.10564v2, section 7; published as *On
Palindromic Periodicities*, CPM 2025, LIPIcs 331, article 11, section 7.

Run `python3 check.py`. The script uses only Python's standard library, integer
arithmetic, and exact rational output from `fractions.Fraction`. Its saved stdout
is `output.txt`; successful execution ends with `ALL_CHECKS_OK`.

The script reproduces all 21 counts in the preprint and the next three terms of
OEIS A374495. It independently recognizes every binary word of lengths 1 through
14 using the literal definition, including either empty palindrome, and compares
these sets with generation by nonempty symmetric roots and cuts. It also checks
four source examples.

The proposed answer is A(2k) ~ 10k 2^k and A(2k+1) ~ 7(2k+1) 2^k. Supporting
checks cover presentation cardinalities through length 24; 34,884 pairs of
presentations through length 16; explicit intersection sets through length 8;
84 cutoff/Bonferroni cases for 4 <= n <= 24 and 0 <= D <= 3; and both exact
presentation-sum identities for 1 <= k <= 100. Pair intersections are counted by
the connected components of the full equality graph, while the proposed bound
retains only the two palindrome matchings. Counts concern distinct words.

These finite checks test source fidelity and the proposed proof's intermediate
identities and inequalities. They do not prove the infinite asymptotic statement.
The asymptotic ratios in the output are exact fractions, with no acceptance
threshold inferred from them. No Lean proof is included.
