# Hanna A384819: exact finite prime-power checks

Run `python3 check.py` with Python 3; no third-party packages or network access
are needed. `output.txt` is its complete deterministic output.

The script independently constructs the bounded correction by the coefficient
recurrence for the formal exponential and by a proper-divisor Euler-ledger
recursion. It compares all 520 terms in Paul D. Hanna's published A384819 table,
and all 41 displayed coefficients of A384820. It checks every prime power at
most 2048, the ghost divisor-sum and prime-power difference identities, an exact
Euler product through degree 96, and a direct rational Taylor exponential through
degree 32. The source-range check covers 117 prime-power indices. Exhaustive
finite interval uniqueness and a deliberately nonintegral correction provide
additional controls. These are finite checks, not a proof of the conjecture.

Sources: [OEIS A384819](https://oeis.org/A384819), revision 14, and its
[b-file](https://oeis.org/A384819/b384819.txt); [OEIS A384820](https://oeis.org/A384820),
revision 12. The literal tables embedded in the script are attributed to Paul D.
Hanna and The OEIS Foundation Inc., under the
[OEIS End-User License Agreement](https://oeis.org/LICENSE).

`check.py` SHA-256: `f36ce92f3986ede4ecae3c75d79e10caf2bd083ddb5e5d7994412da3c92bc86b`.
