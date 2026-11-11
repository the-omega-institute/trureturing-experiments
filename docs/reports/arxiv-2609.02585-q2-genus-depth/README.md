# Exact checks for Pastorek Question 2

Source: Ján Pastorek, *Extremal Asymmetric Depth of Planar Graphs and Hidden Near-Mirror Symmetries of IPR Fullerenes*, arXiv:2609.02585v1, submitted 2 September 2026. Question 2 asks separately about exact attainment of the Heawood bound and square-root growth in orientable genus.

Run with Python 3.12 or newer:

```sh
python3 docs/reports/arxiv-2609.02585-q2-genus-depth/check.py
```

The program is original code using only Python's standard library. It uses integers and exact rational arithmetic throughout. `output.txt` contains its output and ends in `ALL_CHECKS_OK`.

- Fidelity: the Appendix A graph has 20 vertices and 30 edges, is connected and cubic with girth 5, and its stated nonidentity partial automorphism has rank 17. The graph and map are transcribed from the paper, rather than from its software repository.
- Small graphs: all 1,098 labelled simple graphs on 2–5 vertices (771 connected) satisfy the checked local-transposition and averaging bounds, with asymmetric depth computed directly from induced-subgraph isomorphisms.
- Partial injections: all 145,793 nonidentity partial injections on 2–7 vertices satisfy the independent-edge-constraint bound. Domains and images can differ. For all 222 maps on 2–4 vertices, direct enumeration of edge assignments also agrees with the exact constraint-rank probability.
- Genus scalars: the radical comparison holds for genera 1–100,000 by integer squaring. Its two sides are equal at genus 1. For genus at least 3 the general comparison has the positive polynomial certificate `1296 t² + 5040 t + 3616`, with `t = g - 3`; genus 2 is separate.
- Existence budget: the weighted-binomial inequality `33^100 * 32 < 2^10 * 32^100`, the 64 base cases and induction step for `n² ≤ 2^floor(n/64)`, and the 960 residue cases for the exponent estimates support a general union-bound argument for every `n ≥ 65536`. Exact budget values are also checked at five sample orders.

The proposed general counting argument is as follows. A partial injection of rank `r` moving `s` domain vertices changes at least `s(2r-s-2)/2` unordered pairs. Its graph of equality constraints has maximum directed indegree and outdegree one; each nontrivial cycle loses only one independent constraint. Hence it imposes at least `s(2r-s-2)/4` independent fair-edge-bit equations. For omission budgets at most `floor(n/100)`, this is at least `s n/5`. Counting the moved inputs and their images by `n^(2s)` and summing over all omitted domains gives failure probability at most

`2^ceil(n/10) / (2^floor(n/6) - 1)`.

A disconnected random graph has probability at most

`n / (2^floor(n/2) - n)`.

Both terms are below `1/4` throughout the stated range. The first bound follows from the residue estimates and the second from its two parity base cases and monotonicity under `n → n+2`. Thus the main positive result is a general probabilistic argument. The finite tests are supporting checks, not an exhaustive construction of the large graphs or a formal proof of existence.

The genus transfer additionally uses the classical Euler edge bound for minimum orientable genus and the standard rotation-system interpretation of orientable cellular embeddings. Neither topological fact is established by this check program. No Lean proof of the settlement is included.
