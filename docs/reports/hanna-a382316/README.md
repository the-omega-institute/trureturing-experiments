# Hanna A382316 coefficient positivity

Run `python3 check.py` with Python 3. The program uses integers only and has no third-party dependencies. Its output is in `check.out`.

Source: Paul D. Hanna, [OEIS A382316](https://oeis.org/A382316), revision 11 (2025-05-10). The input `b382316.txt` is Hanna's [published table](https://oeis.org/A382316/b382316.txt) of terms 1 through 1030. Source attribution and reuse terms: [OEIS end-user license](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement).

The check reproduces every published term using two different exact recurrences: the positive coupled system and the original denominator-cleared equation. Independent truncated-series multiplication and unit inversion check the quotient bridge, both coupled identities, and the source identity. It also checks integrality and strict positivity on this finite range. This supports source fidelity and indexing; it proves no assertion for all indices.
