# Butterfly gradient variances at half filling and four particles

Haobo Ma and Wenlin Zhang, 18 October 2026.

The note proves an exact first-phase gradient variance formula for the
Walsh-spread magic-block input in one increasing-stride butterfly with
independent uniform phase and rotation angles. At half filling, k=n/2 and
n=4*2^r, all six readouts within modes 0 through 3 have

    V/H ~ (273/224)*(9/8)^r.

Thus the full range including half filling has no uniform finite upper
constant. The sparse uniform comparison remains open. The note also
completes the exact k=4 formulas for these six readouts, their two-sided
constant comparison and their leading relative deficit. Earlier same-thread
k=2 calculations are retained in the coefficient certificate.

## Files

- `note.pdf`: six-page mathematical note with proofs and verified source citations.
- `note.tex`: LaTeX source.
- `verification.zip`: all main verification programs, exact certificates and
  tables, and separate referee implementations with their check receipts.

## Reproduce the verification

Extract `verification.zip` into an empty directory. It contains `verification/`
and `review-work/`, along with this README and the manuscript source.
From the extraction root, run:

```sh
cd verification
PYTHONDONTWRITEBYTECODE=1 python3 code/verify_certificate.py
PYTHONDONTWRITEBYTECODE=1 python3 code/check_round1.py
PYTHONDONTWRITEBYTECODE=1 python3 code/check_haar.py
```

The main checks require Python 3 and its standard library. See
`verification/README.md` for the logical scope and commands that regenerate
every direct integration table. There are 18,918 initial-value/recurrence
identity pairs, six all-depth half-filling transfer certificates, 6,603
full-Walsh RDM entry comparisons, 36 Fourier/recurrence cross-checks, and
120 exact finite-size table entries for pairs 01 and 13 at n=16,32,64,128.

The optional referee layer uses separate implementations and Python + SymPy.
From the extraction root, run:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 review-work/independent_check.py
PYTHONDONTWRITEBYTECODE=1 python3 review-work/independent_transfer.py
PYTHONDONTWRITEBYTECODE=1 python3 review-work/independent_certificate.py
PYTHONDONTWRITEBYTECODE=1 python3 review-work/check_tables.py
```

These independently rebuild Fock RDMs, backwards-path phase/angle integrals,
all sixteen five-state transfers and their low-layer terminal vectors, and
the entire coefficient-certificate verification. They share the delivered
certificate data, not its solver. The finite-table referee recomputes
selected data and checks the remaining recorded fractions and endpoints;
it does not independently reintegrate all 120 large-size entries.

All mathematical computations use exact rational arithmetic. Decimal fields
are for display. These checks are not Lean formalizations. The primary source
is I. Kerenidis, *Scalable Quantum Machine Learning: Trainability, Expressivity
and Efficiency*, arXiv:2607.24014v2 (2026), especially Theorem 9, Claim 1,
Theorem 19 and Appendix A.2: https://arxiv.org/abs/2607.24014v2.

Suggested citation: H. Ma and W. Zhang, *Butterfly gradient variances at half
filling and four particles*, research note, version dated 18 October 2026.
Cite the pinned repository revision containing this directory.
