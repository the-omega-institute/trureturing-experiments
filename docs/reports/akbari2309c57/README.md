# Akbari–Chuharski–Sharan–Slonim, Conjecture 5.7

Run `python3 check.py` from this directory with Python 3.12 or newer.
[check.out](check.out) records the expected output.

The program computes generated subgroups and tests solubility through their
derived series. It reproduces the involution result for all 15 involutions
in the alternating group of degree five (Corollary 5.6(1)). For the products
with cyclic groups of orders 3, 5, 7, and 9, it checks the solubilizer of the
specified element against the full product fiber. It computes radical
membership through normal closures. It also tests the order-four element
in the determinant-one matrix group over the field of five elements.
All comparisons use exact arithmetic.

The odd-order cases in the output sample the proposed infinite family.
The general construction uses the invariance of the solubilizer probability
under quotienting by a soluble normal subgroup.
[Source, arXiv:2309.09104v3](https://arxiv.org/abs/2309.09104v3)
