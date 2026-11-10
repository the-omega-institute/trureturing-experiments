# Exact checks for Conjecture 2.4

Source: Siva Athreya, Satyaki Mukherjee and Soumendu Sundar Mukherjee, *On a generalisation of the coupon collector problem*, arXiv:2304.01145v2, Conjecture 2.4.

Run `python3 docs/reports/2304.01145-conjecture-2.4/check.py` with SymPy installed. The stored output is `check.out` and ends with `ALL_CHECKS_OK`. All evidence uses integers, `fractions.Fraction`, and exact SymPy algebra; no floating point is used.

The checks compare the overlap passage recurrence with rational first-step linear systems in 26 cases: `(r,s)` equal to `(2,1)`, `(3,2)`, or `(3,1)`, with `max(4,r+s) ≤ n ≤ 12`. Every member/outsider replacement is independently enumerated on 1,706 actual r-subsets, verifying the overlap transition counts and zero residuals for the claimed hitting means. Symbolic polynomial leading coefficients are checked for those three pairs and `(4,2)`, and exact normalized passage values are printed.

For five small parameter triples the complete coverage process, with collected-label masks, is solved by rational linear systems and its mean is checked against both record-cover bounds. The logarithmic lower estimate uses an upper enclosure for the logarithm: the rational atanh series plus a geometric remainder bound. Thus the certified rational inequality implies the actual logarithmic inequality. All 62 ordered group compositions with two through six labels, including a time-zero group, are checked by enumerating permutations; the strict-record average agrees exactly with the tied-group formula and satisfies the harmonic and logarithmic bounds.

Fidelity to the source is checked by reproducing the Johnson graph's degree `r(n-r)` for `r=s=2` and `3 ≤ n ≤ 7`. The source model's undefined outsider choice at `n=r` is excluded. These finite checks support consistency and fidelity; they do not prove the asymptotic conjecture.
