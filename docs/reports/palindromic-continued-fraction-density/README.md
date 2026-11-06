# Exact checks for Freeman's printed density convention

Run with Python 3 using only the standard library:

```sh
python3 docs/reports/palindromic-continued-fraction-density/check.py
```

The source is David M. Freeman, *Generalized Palindromic Continued Fractions*, arXiv:1604.06755v2, sections 2, 3, 6 and 7. Section 6 prints the upper limit of successive backward prefix-length ratios. The target is the first sentence of Conjecture 7.1 under that convention.

The check reproduces the source's rational continued-fraction example, its two-digit m-palindromes and their squares, and 192 Fibonacci-derived palindromic prefixes. Independent rational evaluation is compared with integer transfer matrices on all 3,279 words over {1,2,3} of lengths 1 through 7.

For the proposed counterfamily, set W₀=a and Wₖ₊₁=Wₖ^(k+2) b^(k+1) Wₖ^(k+2). Encode a as (mr,r) and b as (ms,s). The program scans all 5,188 positive prefix lengths of the encoded W₄ for each of 24 triples with 2≤m≤5, 1≤r,s≤3 and r≠s. It checks the actual next qualifying prefix after each designated length. For m=2, r=2, s=3 the checked successive ratios are 1/2, 2/3, 3/4 and 4/5. All comparisons use integers or rational numbers; no floating point arithmetic is evidence.

These are finite checks. They do not certify the infinite extension, eventual nonperiodicity, an infinite limsup, or the conjecture's negation. Those require the construction and argument in the preregistration, followed by Lean formalization. In particular, a few ratios exceeding 1/φ do not refute a limsup bound on their own. The proposed infinite argument produces ratios bounded below by (k+1)/(k+2) at arbitrarily large indices of the complete qualifying-prefix list. This concerns the printed limsup convention; the liminf variant is a different question.

`check.out` is the program's complete output. Success ends with `ALL_CHECKS_OK`.
