# Turau's amnesiac-flooding radius conjecture

`check.py` independently implements Algorithm 1 of Volker Turau, *Analysis of
Amnesiac Flooding*, arXiv:2002.10752v4. It uses Python's standard library and
integer arithmetic throughout. Run `python3 check.py`; `output.txt` is its output.
The final line is `ALL_CHECKS_OK`.

Fidelity checks reproduce the complete-graph values stated immediately after
Definition 1, single-source termination at eccentricity on paths, and the stated
single-source termination time on odd cycles. The complete-graph check covers
all 1,009 nonempty source sets for orders 3 through 9.

The proposed counterfamily consists of a path with 5k vertices, with new vertices
a_j adjacent to positions 5j+1 and 5j+2 and b_j adjacent to positions 5j+2 and
5j+3. Its k sources are positions 5j+2. For k=1 through 32, exact breadth-first
search checks radius floor(5k/2), path isometry, and distance-two source coverage.
Two implementations (sets of directed edges and incoming bit masks) agree with
the explicit three-round transmission lists. The k=4 graph has 28 vertices,
35 edges, radius 10, and termination time 3 from its four prescribed sources.
Thus its optimized termination time is at most 3, sufficient to refute
k Flood_k(G) >= Rad(G)+k-1, since 12 < 13.

Exhaustive checks of all k-source sets independently establish Flood_k=3 for
k=1 through 4 (21,903 source sets in total). For general k, the paper's proved
Theorem 4(3) together with the three-round execution would imply Flood_k=3;
the script checks that numerical inference for the inspected values. No finite
experiment proves the infinite family, the general lower bound, or novelty.
The preregistration's mathematical argument supplies the proposed general route;
Lean implementation is outside this report.

The program is original code written for this report. No third-party program or
Lean settlement source was copied or adapted.
