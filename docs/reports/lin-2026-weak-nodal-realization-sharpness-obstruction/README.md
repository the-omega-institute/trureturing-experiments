# Weak nodal realization sharpness obstruction

This report checks the universal realization-sharpness reading of the Section 6 question in Jephian C.-H. Lin, *Congruent quotient matrix, nodal domain, and nodal count*, [arXiv:2611.00789v1](https://arxiv.org/abs/2611.00789v1), for the weak-domain bound of Theorem 5.13.

Run with Python 3 and SymPy:

```sh
python3 check.py
```

The stored output uses Python 3.14.8 and SymPy 1.14.0. All matrix calculations use exact integers or rationals. Square roots in the five-vertex spectrum are exact algebraic expressions. The pseudorandom seed is 261100789.

The fidelity check reproduces Examples 5.3 and 5.10, including the spectrum, eigenvector, and weak-domain number 1. The cover enumerator contracts homogeneous and zero-zero edges, enumerates all permitted connected subgraphs, and tests all five clauses of Definition 5.9. It includes omitted zeros and zero-only singleton members, uses sign unanimity for Definition 5.7, and counts distinct pulled-back domain combinations.

For the proposed family on 2r+1 vertices, the program computes W=2 for r=2,3 (12 and 56 combinations respectively). At r=2 the inertia is (2,2,1) and the exact spectrum is -3, -sqrt(3), 0, sqrt(3), 3, so the weak-domain bound has a gap of 1. The Gram identity, characteristic polynomial, and inertia are checked for r=2 through 7.

The realization check tests 72 rational matrices preserving the complete original sign pattern, including zero diagonals, and 72 variants allowing arbitrary left diagonal entries. It checks eigenvector balance, strict diagonal dominance, full row rank, and inertia (r,r,1). For the relaxed variants it independently constructs an exact congruence to r hyperbolic planes and one zero coordinate. Four random original-pattern cases also undergo complete weak-cover enumeration.

These finite exact checks support the general argument proposed for preregistration; they do not constitute a Lean proof of the universal obstruction. The bound itself is not challenged. No Lean settlement proof is included.
