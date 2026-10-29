# Odd-spacing Jahangir spectrum: exact finite checks

Source: Diego Bravo, Florencia Cubría and Annabella Zapattini,
*From Wheel to Jahangir graphs: spectra and energy under structured edge deletion*,
arXiv:2610.28861v1, Problem 3, https://arxiv.org/abs/2610.28861v1.

Run with Python 3 and SymPy (checked with SymPy 1.14.0):

```sh
python3 check.py > output.txt
```

The graph is constructed as an integer adjacency matrix of a cycle of length
`n=m*q` and one hub joined to residues 0 modulo `m`. Quotient matrices are
computed from neighbour counts, verifying every vertex in every distance cell,
and cross-checked by `A*P=P*Q`.

Fidelity checks reproduce the source's m=2 factorization for q=2..8, and its even
factorization (Proposition 4.5 and the ensuing corollary) for m=4,6,8,10 and q=1..6.

The odd grid exhausts m=3,5,...,15, q=1..8, mq<=80 (50 cases). In every case it
checks cycle-polynomial divisibility, the complete characteristic-polynomial
identity, quotient entries, dimension counts and coprimality of the second
inherited factor with the quotient. For r=1..10 it also checks the quotient
formula symbolically in q and the antisymmetric tridiagonal block. The displayed
m=3 and m=5 formulas are checked symbolically. All evidence uses exact arithmetic;
there are no approximate eigenvalues or floating-point assertions.

`output.txt` is the program's output; success requires exit 0 and final line
`ALL_CHECKS_OK`. These are finite checks supporting a proposed uniform argument,
not a proof for all parameters and not Lean kernel verification.
