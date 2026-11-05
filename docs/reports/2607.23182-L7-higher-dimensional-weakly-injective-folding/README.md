# Higher-dimensional weakly injective folding

Source: Pengzhou Wu, *Beyond ICA: Identifiability by Symmetry Breaking*,
arXiv:2607.23182v1, Appendix L.7, with Example A in Appendix L.3 as the fidelity case.

Run with Python 3, using only the standard library:

```sh
python3 docs/reports/2607.23182-L7-higher-dimensional-weakly-injective-folding/check.py
```

The program uses integers and `fractions.Fraction` throughout. It reproduces
Example A's formula, continuity and chamber fibers. It checks 18 rational maps in
dimensions 2 and 3, including explicit matrices with determinants (1,-1,1) and
(3,-4,12). Exact matrix identities establish boundary continuity and the inverse
first-coordinate formulas for every target; 540 rational targets also exercise
all three chambers and both boundary hyperplanes. The chamber fibers are 1,3,1
and the directional active-set differences across successive boundaries are
(0,2) and (2,0), so neither boundary is simple.

The continuity constraints on arbitrary affine triples have rank 2d²; the full
space of triples on this fixed partition has dimension d²+3d (10 and 18 in the
checked dimensions). Strict alternating determinant signs specify the proposed
open family in that space. These checks support the general mathematical
argument; they do not certify the arbitrary-dimensional existence assertion or
any probability or genericity claim. No Lean proof accompanies this report.

`output.txt` is the program's output; its last line is `ALL_CHECKS_OK`.
