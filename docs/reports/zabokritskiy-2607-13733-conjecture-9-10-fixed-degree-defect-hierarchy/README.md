# Fixed-degree defect hierarchy: exact finite checks

Source: Aryeh Lev Zabokritskiy (Yohananov), *Independent Sets in Multiset Profile Graphs via Weighted Local Covers*, arXiv:2607.13733v3, Conjecture 9.10 and Eq. (69).

Run with Python 3.12 or newer, with assertions enabled:

    python3 check.py

The standard-library program uses only integers and `fractions.Fraction` for evidence. It independently builds the source's triangular exact-coverage recursion and the reciprocal-fixed-point average using the derangement distribution. Direct permutation enumeration checks the distribution through size seven. It reproduces the four degree-seven source weights; compares every anchor partition for degrees through fourteen; checks the marked-fixed-point identity for every vertex partition through degree sixteen; and verifies every admissible truncated level through degree sixteen. The last check exhibits the failing degree immediately below the proposed threshold in eight small cases, for triage only.

`output.txt` is the complete deterministic output. The final line is `ALL_CHECKS_OK`. These are finite diagnostic checks, not a proof of the conjecture or its negation. No Lean proof or kernel verification is claimed. The source's displayed asymptotic consequence is outside the settling claim.
