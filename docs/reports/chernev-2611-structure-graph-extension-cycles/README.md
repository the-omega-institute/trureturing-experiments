# Structure-graph cycles over finite extension fields

Source: A. M. Chernev and N. A. Chernev, *On the Connected Components of Structure Graphs of Associative Rings*, arXiv:2611.07521v1, Hypothesis 6.1, Definition 3.1, Example 5.3 and Section 8. The map is Φ(a,b)=(a+b,ab).

Run `python3 check.py`. The standard-library program uses integer coefficient arithmetic in explicit polynomial quotient fields. Every modulus is checked for irreducibility by monic trial division; every multiplication is compared between shift-and-reduce and coefficient convolution followed by polynomial long division. No source program was copied.

`output.txt` contains the exact output. Fidelity checks reproduce the source's Z₅ cycle and the component inequalities for GF(9) and GF(49). The checks also test Φ²(−1,x)=(−1,x−x²), eight iterations of that identity at every element, nonzero periodic-point lifting, odd-degree prime-field-complement invariance, and the odd-characteristic four-cycle wherever applicable. GF(8), GF(9), GF(16), GF(25) and GF(32) are included, along with boundary fields and additional cases. The GF(16) eight-cycle is checked arrow by arrow. A separate integer check tests the characteristic-two degree split through degree 256.

These are finite checks supporting a proposed general argument. They do not establish the universally quantified hypothesis, which requires a subsequent Lean proof. The GF(4) exception has only loops. The companion repository's README and MIT license were read; this program is an independent implementation.
