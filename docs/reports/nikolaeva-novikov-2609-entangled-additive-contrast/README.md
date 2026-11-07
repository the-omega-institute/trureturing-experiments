# Entangled density matrices with additive logarithmic contrast

`check.py` checks the noncommuting entangled family proposed for Nikolaeva–Novikov, arXiv:2609.06997v2, Remark 5.25 and Section 6.3 (third open-question bullet). It uses exact SymPy algebra throughout.

Run `python3 check.py` with Python 3 and SymPy. The recorded output is `output.txt`; success ends with `ALL_CHECKS_OK` and exits 0.

The fidelity checks reproduce the source's Theorem 4.4 matrices, their singular extrema, and the partial-transpose matrix and spectrum in Section 5.9. The family checks use k=3+t with t>0, two explicit orthonormal eigenbases, trace normalization, the common extremal eigenvectors, two negative partial-transpose expectations, and a strictly positive commutator coefficient. Exact principal-minor certificates with nonnegative rational polynomial coefficients bound every middle-block singular value between 1 and k² for all t>0. Thus the product condition number is k², using singular values even though the product is not Hermitian.

There are 37 exact checks. These checks support the proposed analytic argument; they are not a Lean proof. The family is a uniform argument over all real k>3, rather than a finite enumeration of positive instances. No claim is made for k≤3 or for optimality of this threshold.
