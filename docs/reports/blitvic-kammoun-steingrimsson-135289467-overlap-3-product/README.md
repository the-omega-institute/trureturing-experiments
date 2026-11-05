# Pattern 135289467 with overlap 3

Exact finite checks for Natasha Blitvić, Mohammed Slim Kammoun and Einar
Steingrímsson, *Permutations, Pattern Covers, Posets and Positivity: Part I:
Patterns and Posets*, arXiv:2610.29846v1, Section 7 item 10.
Source: https://arxiv.org/abs/2610.29846v1.

Run with Python 3 and no external dependencies:

```sh
python3 check.py > output.txt
```

The program builds the original position-comparison poset from every pair of
entries in each prescribed window. An integer ideal dynamic program counts its
linear extensions without assuming the proposed recurrence. Additional pattern
occurrences are allowed. It checks the target for every j from 1 through 9,
the explicit seven-element head and six-element recursive gluing for every j
from 2 through 12, and the gluing count on all 76 forward-edge DAG presentations
on at most four vertices. Distinct DAG presentations can represent the same
poset; 76 is not a count of poset isomorphism classes.

Fidelity checks reproduce the source's Proposition 5.25 for pattern 125634 and
overlap 2 through j=5 and its Remark 2.7 example 192837465, whose extra occurrence
is permitted. Direct enumeration of all 9! permutations checks the target at
j=1. Polynomial coefficients and all counts use exact integers; no floating
point or numerical approximation occurs.

The observed target counts are
1, 30, 3630, 1049070, 554958030, 466719703230, 571731636456750,
961080880883796750, 2123027665872307020750.
They agree with A_1=1, A_2=30 and A_j=(6j-7)^2 A_(j-1) for j>=3.
The output ends in `ALL_CHECKS_OK`. These are finite checks supporting a
preregistration; they do not prove the universal claim or supply Lean evidence.
