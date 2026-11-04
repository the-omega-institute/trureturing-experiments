# Exact checks for Ye's dense-phase consensus-gap question

Source: Hao Ye, *Self-Attention Dynamics with Rotary Position Embeddings: Twisted States and Explicit Consensus Rates on the Sphere*, arXiv:2607.24502v1, Sections 5 and 9 and Appendix B.1.

Run with Python 3 and SymPy:

```sh
python3 docs/reports/ye-rope-dense-phase-gap/check.py
```

`output.txt` is the program's output. The program uses integers, rational numbers and symbolic exact algebra; it performs no floating-point evaluations.

The fidelity checks recover the complete spectrum and twelve zero eigenvalues of the source's Figure 1(a) ring (n=16, m=4), symbolically for x=exp(beta), including beta=2. They also reproduce Theorem 5.3(ii)'s exact L=2 gap. The route checks verify phase quotient labels for every signed resonance with 2<=n<=32, the proposed cooperative Fourier-difference system for every 4<=L<=64, and the quadrature bound on monomials of degrees 1 through 16 on meshes of 1 through 32 cells.

The symbolic difference-system check includes both parity boundary formulas and nonnegativity after adding one half of the identity. This checks the algebra needed for the proposed all-temperature ordering of Fourier modes. The checks support the proposed general proof; they do not establish the unbounded joint asymptotic. No proof of the claim or its negation is supplied in Lean.
