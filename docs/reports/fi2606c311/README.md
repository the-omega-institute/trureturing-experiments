# Frank–Ivković excess checks

Source: Michael Frank and Stefan Ivković, *Near-Riesz and pseudo-Riesz bases of Hilbert C*-modules and semi-C*-Fredholm operators*, arXiv:2606.16597v3, Definitions 2.3–2.4 and 3.6–3.7, Example 3.3 and Conjecture 3.11. The version was submitted on 4 November 2026; its TeX title date is “November 3rd, 2026, version 5.3”.

Run `python3 check.py` with Python 3.12 or later and SymPy (checked with Python 3.14.8 and SymPy 1.14.0). The stored standard output is `check.out`. Successful execution exits 0 and ends with `ALL_CHECKS_OK`.

The fidelity check reproduces Example 3.3 in the coefficient algebra C²: the frame (p₁,p₂,p₂) has frame operator diag(1,2), canonical coefficients (p₁,p₂/2,p₂/2), and carrier-corrected excess p₂. Its Riesz subframe (p₁,p₂) has excess zero.

For the proposed counterexample, the program checks exact symbolic identities for the frame (t,1), synthesis kernel, kernel projection, canonical coefficients, carrier approximants, and their errors. It also checks 65 rational points with 64 approximants each, 128 values of tₘ=1/m, and the corresponding identities for an arbitrary complex parameter f. All arithmetic is exact; no floating-point computation is used.

The analytic argument still to be formalized concerns the actual universal bidual of C([−1,1],ℂ), the weak-star carrier limit, normal evaluation extensions, and failure of continuity of the excess readings at zero. The program does not certify those steps or prove the original conjecture or its negation. The proposed counterexample targets membership in the coefficient algebra; its excess remains a bounded positive projection. The kernel-projection diagonal sum is 1, while the carrier-corrected excess reading at zero is 0; this is an erratum candidate for Proposition 3.9(ii).
