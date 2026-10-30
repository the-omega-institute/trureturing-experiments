# Source and licensing

Wang–Dong, arXiv:2610.39833v1, Conjecture 6:
https://arxiv.org/abs/2610.39833v1 . For every m >= 4 and admissible n,r with
s=n-m and 1 <= r < s <= n, the binary and arbitrary finite alphabet maxima
both equal min(n-r,3m-1). Automata are total, finite, strongly connected and
synchronizing; S=Qv is an actual reachable image, and A is a nonempty proper
subset of S. Each original letter has unit cost. Success is avoidance
(A not subset Sw) OR compression (|Sw| < s).

The source license is arXiv nonexclusive, as supplied by prior documentary
qualification; it was not independently requalified here. This repository
retains only citation and limited mathematical reference, with no bulk paper
text. No paper rights are inferred from the code license.

The independently written programs are Apache-2.0 under the repository LICENSE.
Z3 / z3-solver==4.15.4.0 is an unvendored MIT-licensed dependency:
https://github.com/Z3Prover/z3 . No solver binaries or dependency sources are
included. `upfront-cuts/encoding.py` derives from `readiness/readiness.py`;
its original source SHA256 is
47b0299e9067063027ece241d33120513f7a26108ee4f60c389aaa44dfc996e1.
