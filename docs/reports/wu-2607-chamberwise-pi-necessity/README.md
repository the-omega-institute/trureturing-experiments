# A chamber-wise PI obstruction to map identifiability

Source: Pengzhou Wu, *Beyond ICA: Identifiability by Symmetry Breaking*,
arXiv:2607.23182v1, Section 5, Theorem 8 and Appendix E, Proposition 18.

Run `python3 check.py` with Python 3 and SymPy. The output is in `check.out`.
The checker uses integer and rational arithmetic and symbolic polynomial algebra.
It reproduces Proposition 18's global parameter collision and chamber-wise
injectivity, then checks the proposed two-dimensional reflection construction.

All 64 rectangular cells are covered, including the unbounded cells. Each sign
matrix sends an entire input cell onto a computed output cell. The two maps have
four distinct nonempty affine rules. Every output cell has chamber-wise PI and
matching weighted Gaussian parameter multisets. The program checks the four
component reflection identities, the localized-symmetry composition identities,
all branch USB witnesses, positive definiteness, and the covariance/mean algebra
used in mixture symmetry triviality. Both maps have two preimages throughout the
regular part of the upper square. Every cell and the final checks are printed.

The proof obligation of finite Gaussian-mixture identifiability is separate from
the verified finite algebra. Equality of kernel sums is also separate from the
measure-theoretic pushforward proof. No Lean proof is provided by this check.
The proposed settlement concerns **MID** (Theorem 8). The two latent laws are
identical, so **LID** (Theorem 5) holds via the identity affine map.

The same componentwise identities hold for weights lambda and 1-lambda. The
stated family uses 0 < lambda < 1 and lambda != 1/2, so distinct weights fix the
component matching in the MST argument. The main planned proof uses localized
measure-preserving symmetries; finite cell computations are helper checks.
