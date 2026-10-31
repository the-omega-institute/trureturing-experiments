# Corrected Mathar–Yuen recurrence for OEIS A200770

Source: [OEIS A200770, revision 14](https://oeis.org/A200770/internal), with
Hardin's array definition and Mathar's FORMULA conjecture corrected by Yuen.
The program embeds the 41 printed terms and all seven length-4 examples.

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The stored stdout is
`output.txt`; successful execution exits 0 and ends with `ALL_CHECKS_OK`.
All evidence uses integers, rational fractions or symbolic integer polynomials.

The check enumerates all 265,720 ternary words of length at most 11, comparing
literal nondecreasing rational prefix averages with the proposed normal form.
It checks the fixed-content count for all 364 content triples in that range,
the balanced-content formulas through index 500, the recurrence against every
applicable printed term and through index 500 using the formulas, and three
symbolic identities. The output supplies all 15 numerator polynomials and their
zero coefficient sums. This is evidence for a proposed all-index argument;
it is not a Lean proof of the conjecture or of the word-classification lemma.

Proposed mechanism: with r zeros, s ones and t twos, zeros are initial, and a
one can occur only before the (r+1)-st two. Thus the word has the unique form
0^r v 2^(t-min(r,t)), where v has s ones and min(r,t) twos. Choosing the one
positions gives binomial(s+min(r,t),s). Summing balanced contents produces the
three residue formulas consumed by the symbolic checks.
