# Exact checks for logarithmic-lattice classification

Source: Ciro S. Campolina and Alexei A. Mailybaev, *Fluid dynamics on logarithmic lattices*, arXiv:2005.14027v1, Theorem 1, Remark 2, Lemma 3 and Table 1. The check concerns the conjectured removal of the spacing and exponent cutoffs. The publication DOI is https://doi.org/10.1088/1361-6544/abef73.

Run with Python 3 and SymPy:

```sh
python3 check.py > check.out
```

The program uses integers, rational arithmetic and exact polynomial factorization and root counts. It reproduces the source's three dyadic and twelve ordered plastic-number unity triads, as well as its exact cutoff inequality at 21/20. It factors every polynomial X^n-X^m-1 with 1 <= n <= 60 and 0 <= m < n, reconstructs every factorization, checks the claimed cyclotomic factor and irreducible quotient, and groups identical noncyclotomic factors. The only shared positive-root pairs in this finite range are (3g,g) and (5g,4g), for 1 <= g <= 12. The 60 cases m=0 are included.

The component checks concern finite cyclic projections of the sign/exponent triad graph. They include all three sign patterns coming from a normalized positive relation and compare breadth-first search against the predicted gcd residue classes. They are checks of a combinatorial model; the bridge from real lattice triads and the statement on the infinite lattice still require a general proof.

`check.out` is the output of the program in this directory. Its final line is `ALL_CHECKS_OK`. These finite checks neither prove the unrestricted conjecture nor certify any Lean theorem. No floating point calculation is used as evidence.

The general proposed proof compares irreducible polynomial factors using the Ljunggren–Tverberg factorization, then compares the first nonconstant coefficient, and finally uses integer gcd connectivity. The finite computation is only an independent consistency check of that argument.
