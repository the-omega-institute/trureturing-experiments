# Sain–Paul–Debnath, Conjecture 2.11

Run `python3 docs/reports/sain1407c211/check.py` from the experiments repository root. The captured output uses Python 3.14.8 and SymPy 1.14.0.

The fidelity check reproduces the strongly orthonormal cube basis in Example 2.9 and exhibits a vector in the complementary span which disproves zero best coapproximation. The family check uses exact rational and radical arithmetic for conjugate planar norms and their rotational lifts. It checks values and gradients at seams, agreement of conjugate meridian directions, all entries of the supporting-functional pairing matrix, and the parallelogram defect. The sampled cases include the Euclidean control, both non-Euclidean sector orderings, the rotation axis, opposite signs of the axial coordinate, and a horizontal frame with non-coordinate axes.

These checks cover the listed finite sample. The analytic argument for the continuum of exponents, dimensions, and prescribed unit vectors belongs to the preregistration draft.

[Source paper, arXiv:1407.5016v3](https://arxiv.org/abs/1407.5016v3)
