# Kurkov A107667 weighted-ballot checks

Run `python3 check.py` with Python 3. The program uses only the standard library and exact integers and rational numbers. `output.txt` contains the resulting output.

The source is [OEIS A107667](https://oeis.org/A107667), approved revision 25 (7 June 2026); the [history](https://oeis.org/history?seq=A107667) places Kurkov’s conjecture in revision 23 (2 June 2026). Paul D. Hanna originated the entry. Source data are attributed to that entry under the [OEIS End-User License Agreement](https://oeis.org/LICENSE).

Fidelity checks reproduce all eight explicitly printed triangle rows, all seven explicitly printed square rows, and the source matrix diagonalization through row 15. The check distinguishes the earlier conjugating matrix from Kurkov’s three-argument table, although both are named P in the source.

The native shifted-square recurrence, Kurkov’s printed P recurrence, and forward weighted-path dynamic programming agree on all 861 entries with `0 <= b <= a <= 40`. The program also checks 12,341 P entries, 820 first-return convolutions and 54 direct path enumerations. These are finite computations; they do not prove the all-index conjecture and are not a Lean kernel check.
