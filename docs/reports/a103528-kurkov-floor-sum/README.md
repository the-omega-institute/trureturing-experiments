# OEIS A103528 floor-sum difference

This report checks Mikhail Kurkov’s 1 June 2022 FORMULA conjecture in OEIS A103528. The program uses Python integers only: `A(n)` is the defining congruence sum, `B(n)` is the finite floor sum with the exact base-two logarithm cutoff, and every case n = 2 through 4096 is checked against `2*A(n) = B(n) - B(n-1) - 1`. The prefix in the OEIS entry is checked as a fidelity result.

Run `python3 check.py`; the final line is `ALL_CHECKS_OK` on success. The arithmetic uses no floating point.
