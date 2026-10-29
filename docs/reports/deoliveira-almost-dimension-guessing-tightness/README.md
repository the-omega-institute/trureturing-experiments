# Almost-dimension guessing-bound tightness: exact finite checks

Run `python3 check.py` with Python 3 and SymPy 1.14.0. The stored output is `check.out`.

The source is A. de Oliveira Junior et al., *Quantum prepare-and-measure scenarios: concepts and methods*, arXiv:2610.20725v1, Eq. (221), and J. Pauwels, S. Pironio and A. Tavakoli, *Information capacity of quantum communication under natural physical assumptions*, Quantum 9, 1637, arXiv:2405.07231v3, Eqs. (6)–(8), DOI 10.22331/q-2025-02-18-1637.

The fidelity check reproduces the source's vacuum-component envelope from the equiangular Gram spectrum and pretty-good-measurement expression, including the reduction of Eq. (8) at d=1. This means fidelity to the source, not tightness for arbitrary fixed fidelity targets.

The construction check uses a primitive nth root of unity in the exact quotient by its cyclotomic polynomial. It verifies both Fourier orthogonality and completeness of the rank-one measurement, then checks normalization, every input's common-subspace weight and the exact success formula for n=2 through 9, every 1<=d<n, five leakage values up to the threshold and two values on the plateau. Rank-one outer products are positive semidefinite by their factorization. A separate four-dimensional matrix calculation gives projected weight 9/10 and success 4/5 at d=2 and epsilon=1/10.

All arithmetic is exact; no floating point or SDP solver supplies evidence. These finite checks support a preregistration outline and do not establish a universal theorem or a Lean proof. The program is independently written for this experiment; no third-party source code is copied.
