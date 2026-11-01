# Feldman Conjecture 9.2: exact preregistration checks

Source: David Victor Feldman, *The Missing Zigzag: Cycles of Semitone Trichords and a Conservation Law in Equal Temperament*, arXiv:2609.26114v1, submitted 21 August 2026, Section 5 equation (2) and Conjectures 9.2–9.3. The paper and documentation are CC BY 4.0. The checker is an independent implementation of the displayed mathematical definitions; no author code is copied.

Run `python3 docs/reports/feldman-zigzag-conjecture-nine-two/check.py` from the repository root. It uses Python 3's standard library, integers and rational numbers. The stored output is `check.out`.

The checker counts actual labelled choices by matching opposite full imbalance vectors for t = 2, 3, 4, 5 (n = 6, 9, 12, 15); separate whole-choice enumeration cross-checks t = 2, 3. It preserves coincident endpoint labels, permits only forms II–V at class 2, and imposes balance without closure or connectivity. These counts equal 4, 12, 48, 216.

It reproduces the paper's complete t = 2 through 12 table using the repository's existing even/odd closed forms, checks both exact ratio and error identities for r = 1 through 200, and checks the uniform error bound through t = 400. The owners of those all-parameter closed forms are `D5/S3/Combinatorics/Zigzag/FeldmanConjectureNineThree.lean` in `the-omega-institute/trureturing`. The table reproduction beyond t = 5 uses those closed forms rather than independent enumeration.

All arithmetic used as evidence is exact. Finite checks support the preregistered general argument and do not prove the unbounded limit. A passing run ends with `ALL_CHECKS_OK`.
