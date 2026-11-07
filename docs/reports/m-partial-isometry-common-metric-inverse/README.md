# A prescribed inverse with no common positive metric

The question is Aouichaoui–Buchała–Garcia, *On m-partial isometries: spectra,
weighted shifts, and similarity*, [arXiv:2606.12865v1](https://arxiv.org/abs/2606.12865v1),
Section 6, Question 6. The adjoints in S = 2T* − (T*)²T are fixed in the original
inner product. Both Penrose symmetry conditions must hold in the same new inner
product.

Run with Python 3 and SymPy:

```sh
python3 check.py
```

The saved output is `output.txt`. Success ends with `ALL_CHECKS_OK` and exit 0.
The program uses rational and symbolic exact arithmetic throughout.

For fidelity, it applies the source's Lemma 2.8 to the source's own example before
Question 6, reproducing its displayed TS matrix, TST = T, STS = S, and
ST = I − beta_2(T).

For every real r ≥ 2, let a = sqrt(2 − 1/r), b = (r − 1)/sqrt(r), and let T have
first column (a,b) and zero second column. The checks establish beta_2(T) =
diag(0,1), strictness, and the reflexive inverse identities. They also establish
ST = diag(1,0), (TS)_12 = 2ab, and (TS)_21 = ab(2 − r). Any positive Hermitian
Gram matrix making ST selfadjoint must be diag(p,q), with p,q > 0. Making TS
selfadjoint would require 2p = (2 − r)q. This is impossible throughout r ≥ 2.
The code substitutes r = 2+t with t ≥ 0, and checks that the symmetry residual is
ab(2p+tq) > 0. The endpoint r = 2 is included.

The exact rational witness r = 25 is T = [[7/5,0],[24/5,0]]. Its identities are
also checked independently with Python's fractions.Fraction. Each idempotent is
individually orthogonalizable: Q = ST is already selfadjoint; for P = TS,
W = P*P + (I−P)*(I−P) is positive definite and WP = P*W.

These calculations support a preregistered refutation outline. They provide no
Lean proof, kernel verification, or exhaustive prior-art claim.
