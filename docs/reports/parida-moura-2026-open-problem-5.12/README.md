# Parida–Moura Open Problem 5.12

Exact checks for the matching construction proposed for arXiv:2605.12634v3,
“Cover-free families on graphs”, Section 5, Open Problem 5.12.

Run with Python 3.12 or newer, with no additional packages:

```sh
python3 check.py
```

`check-output.txt` is the corresponding output. Success ends in `ALL_CHECKS_OK`.
The program uses integers and finite sets throughout.

The fidelity check uses the eight blocks printed in Proposition 5.11, verifies
all 24 edge/outside-column relations and both directions of edge separation,
and checks the Sperner capacity bounds giving t(E₈)=5 and t_e(E₈)=4.
The repeated-column ECFF construction on four rows supplies the latter upper bound.

For each d from 3 through 13, let r=floor(d/2) and m=binomial(d,r)+1.
The program checks the proposed matching-CFF on d+2 points: one edge
(S∪{p}, S∪{q}) for each r-subset S of a d-point core, and the additional
edge ({0,…,r}, {0,…,r−1}∪{r+1}). It checks 8,203,100 outside-column relations,
the common block size, and both edge-separation directions. A separate
implementation using sets checks d=3,…,7 (3,688 outside-column relations).
A deliberately broken construction is rejected as a negative control.

For each tested d, exact middle-binomial capacity comparisons give the
Sperner lower bounds d+2 and d+1. A labelled ECFF with identical columns
within each edge on d+1 rows is checked explicitly. These lower bounds use
Sperner's theorem: all CFF columns form an antichain, and all ECFF edge
unions form an antichain. The computations corroborate the uniform argument;
they do not prove an assertion over all d and are not Lean proofs.

The programs implement the mathematics independently in this research seat;
they do not copy third-party program source.

The standalone `cff_verify.py` entry point reuses the bit-mask checker and emits
`cff-verification.json`; run it with `python3 -B cff_verify.py`.
