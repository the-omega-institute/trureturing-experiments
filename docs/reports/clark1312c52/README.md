# Clark–Elhamdadi–Saito–Yeatman Conjecture 5.2

Run `python3 check.py` with Python 3.12 or newer. The program uses integer
arithmetic and the Python standard library. Its output is in `output.txt`.

The trefoil calculation reproduces the source example with 12 nontrivial
colorings in the four-element Alexander quandle. For the proposed family
on ℤ/p² ⊕ ℤ/p, the program checks the target automorphisms, connectedness
and the Seifert determinant. It computes the original two off-diagonal
kernel sizes for each action by elimination modulo p² and modulo p.
At p = 5 it cross-checks those sizes by counting all triples in the
original coefficient group through a table of the third summand.
These checks cover the finite algebra. The coloring interpretation and
prime-knot realization use Liu’s Lemma 5.1 and Friedl’s Theorem 1.1.

- [Clark–Elhamdadi–Saito–Yeatman, arXiv:1312.3307v2](https://arxiv.org/abs/1312.3307v2)
- [Liu, arXiv:2610.34831v1](https://arxiv.org/abs/2610.34831v1)
- [Friedl, arXiv:0704.3390v1](https://arxiv.org/abs/0704.3390v1)
