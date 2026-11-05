# Klech maximum spectral radius with m+2 arcs

Source: Rostislav Klech, *Generating Functions and the Minimum Spectral Radius in Strongly Connected Digraphs with m+2 Edges*, arXiv:2609.18367v1, Conjecture 5.15 and Theorem 5.14.

Run with Python 3 and SymPy 1.14.0:

```sh
python3 check.py
```

The program enumerates every labelled strongly connected 0/1 digraph with m+2 arcs for m=3,4. It compares spectral radii using integer characteristic polynomials and exact algebraic root isolation. Theorem 5.14 supplies the source fidelity check: the smallest spectral radius is the largest real root of X^m-2X-1. The candidate checks verify both maximum values and all labelled equality cases in those sizes (the loopless clause starts at m=4). It independently computes the HH determinant formula for all five listed tuple representations where applicable, for 3<=m<=12, and checks one endpoint identity.

`output.txt` contains the successful output. No floating point is used as evidence. These are finite diagnostics, not a proof of either unbounded clause. The proposed settlement requires the general branch-path argument.
