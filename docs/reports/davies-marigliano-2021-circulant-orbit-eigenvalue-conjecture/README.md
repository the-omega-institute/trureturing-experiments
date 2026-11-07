# Circulant graph pair orbits and adjacency eigenvalues

Source: I. Davies and O. Marigliano, *Coloured Graphical Models and their
Symmetries*, Le Matematiche 76(2), 501–515 (2021), Conjecture 4.4;
[arXiv:2012.01905v5](https://arxiv.org/abs/2012.01905v5).
Unordered vertex pairs include diagonals, as in the source's Corollary 3.3.

Run with Python 3, SymPy 1.14.0 and NetworkX 3.7:

```sh
python3 check.py > output.txt
```

`check.py` uses integers and exact polynomial arithmetic. Its source-fidelity
check reproduces Theorem 3.4(1) for cycles of orders 3–10, including the exact
six-cycle spectrum in Example 3.5. Full graph-automorphism enumeration, without
assuming a product decomposition, gives eight orbits on unordered vertex pairs
for `K₂ □ K₃ □ K₅`; the exact characteristic polynomial gives seven distinct
eigenvalues. The enumerator returns 1,440 automorphisms. The diagonal orbit has
30 pairs and is included in the count.

Direct adjacency characteristic-polynomial checks also cover `K₂ □ K_b □ K_(b+2)`
for `b = 5, 7`, together with connectivity, degrees, common-neighbor counts on
edges, and a cyclic automorphism visiting every vertex. Separate finite
subset-sum checks cover every odd `b` from 3 through 101 and five growing lists
of coprime factor sizes. The output labels counts obtained from a proposed
orbit formula separately from counts obtained by exhaustive automorphism
search. These finite checks support the accompanying infinite-family argument;
they do not replace its proof or a future Lean kernel check.

`output.txt` contains the exact results and ends with `ALL_CHECKS_OK`.
