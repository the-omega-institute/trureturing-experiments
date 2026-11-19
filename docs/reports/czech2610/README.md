# Integer-index monotonicity for qubit stabilizer states

Run `python3 check.py` with Python 3.12 or later. The script uses only the
standard library. The saved output is in `check.out`.

The fidelity check enumerates the Bell and three-party GHZ supports on the
replica torus and recovers their entropy formulas in Table 2 of Czech, Feng,
Harper, Ju, and Zhao, arXiv:2610.26814v1. For the four-party witness, the script
finds parity checks directly from the support words and ranks the full binary
constraint matrix on ket and shifted bra labels. It checks all nontrivial
parameters in the fields of four and eight elements, both replica indices for
each, normalization, and the six two-party support projections. Outputs give
exact contractions and entropy coefficients in units of log 2.

The checks are finite evidence for the preregistration. The draft gives the
argument for the family of larger fields.

[Source, Definition 2.1 and Table 2](https://arxiv.org/html/2610.26814v1)
