# Bowden–Tymoczko Question 6.3

`check.py` performs exact integer enumeration for the minimum-rank part of Question 6.3 from arXiv:1501.02027v1. A divisor `d` of `m` represents the ideal `(d)` in `Z/mZ`. For every triangle, `K4`, path on 2–5 vertices, and cycle on 3–5 vertices, and each `m` in `{4, 6, 12, 30}`, it enumerates all nonzero and all proper nonzero ideal labelings, computes the CRT component rank, and checks the partition-envelope formula. It also checks the source’s rank-one `C3` example with labels `6, 15, 10` modulo `30` by enumerating the spline group and its p-torsion.

All arithmetic is exact integer arithmetic. The output is stored in `check.out`; its final line is `ALL_CHECKS_OK`.
