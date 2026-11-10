# Disjoint critical diameter configurations

Source: Gyivan Lopez-Campos, *On Diameter Graphs Problems: A proof of Brass's Conjecture*, [arXiv:2611.05631v1](https://arxiv.org/abs/2611.05631v1), Section 5, Question 4 (printed page 14). Vertex-criticality is the definition in Section 2.1. The question is stated for general dimension d.

Run `python3 check.py`. The program uses only Python's standard library and exact integer or rational arithmetic. Its successful output is in [output.txt](output.txt), ending in `ALL_CHECKS_OK`.

The source fidelity check reproduces Proposition 1 (existence of a vertex-k-critical subgraph) on all 1,100 labelled simple graphs with at most five vertices. On those graphs it also checks the critical-packing implication: chromatic number at least t*k permits t disjoint induced k-critical subgraphs.

The geometry check enumerates every 2p-subset of [4p] containing the first element, for p = 1, 2, 3. With sign vectors x_A and coordinates z_A(i,j) = x_A(i)x_A(j)/(4p), it checks distinctness, squared diameter one and the exact distance identity

\[
\|z_A-z_B\|^2=1-(|A\cap B|-p)^2/p^2.
\]

For the prime cases p = 2, 3 it also checks every tested separating-polynomial evaluation and its squarefree multilinear reduction over the finite field. The p = 1 geometry case makes no finite-field assertion.

At p = 37, d = binom(148,2) = 10878, the program verifies the strict rational inequality 3(27/16)^37/(2*149) > 2*d. It checks the binomial estimates at every prime from 37 through 100. It also computes and checks the positive coefficients of the polynomial identity supplying the proposed induction step for all p >= 37.

A reuse-interface check evaluates the different classical configuration in `Borsuk.borsuk_theorem` from [mo271/FormalBook](https://github.com/mo271/FormalBook/blob/7df264a4d5e57f9036acdf85f7a2702b3ad82a06/FormalBook/Chapter_18.lean). At p = 37 its source dimension is 10585, which embeds isometrically into dimension 10878, and the ceiling of its partition lower bound is 3137001950. This Python check does not compile or certify that external Lean source.

These are finite diagnostics and exact numerical interface checks. They do not construct the enormous p = 37 configuration explicitly, certify its chromatic number by enumeration, prove the infinite family, settle Question 4 in Lean, or answer the motivating dimension-4 question. The proposed mathematical route combines the classical Kahn–Kalai / Frankl–Wilson input with vertex-critical packing. No third-party program or Lean code is copied into this report.
