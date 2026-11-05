# Turán maximum for the augmented Sombor index

Source: Chunlei Xu, Kinkar Chandra Das and Jayanta Bera, *Structural Properties and Applications of the Augmented Sombor Index*, arXiv:2606.26509v2, §4, the extremal sentence preceding Conjecture 4.1. The target is ASO(G) ≤ ASO(T_n(k)), with equality exactly for graphs isomorphic to T_n(k), on simple k-partite graphs of order n ≥ 3 without a two-vertex path component.

Run with Python 3 and SymPy:

```sh
python3 check.py
```

The program uses integers, rational numbers and symbolic algebra. It normalizes square roots into squarefree radicals and certifies comparisons with rational intervals computed by integer square roots. It never uses floating point as evidence. The captured standard output is in [output.txt](output.txt).

The fidelity check reconstructs the source's complete-bipartite formula from graph degrees and checks its balanced maximum and equality case for every split with 3 ≤ n ≤ 30. The structural conjecture and its equality case are checked for all positive part-size partitions with 3 ≤ n ≤ 20, every balancing transfer among those partitions, and every labeled simple graph with 3 ≤ n ≤ 5 in the source's domain. The graph checks include isolates and disconnected graphs.

The symbolic checks verify the aggregate-curvature identity and the polynomial used in the general proof outline. With R = (4 − 2t)((1+t)⁻² + m⁻²), m ≥ 5, the useful inequality is

(1+t)²((1+2t) − 4t(R−1)) ≥ 1 − 216t/25 + 501t²/25 + 6t³ + 8t⁴/25 > 0.

The quadratic part has discriminant −3444/625. Squaring R−1 gives a false estimate: at m = 5 and t = 1/5, (1+2t) − 4t(R−1)² = −59546/78125. The program explicitly checks that distinction. It also verifies the exact lower-bound polynomial and R(3/4) ≤ 449/490.

These results support a preregistration draft. Finite checks do not prove the universal conjecture. No Lean proof is included. The printed bound in Conjecture 4.1 has two denominator discrepancies relative to ASO(T_n(k)); this experiment targets the structural extremum rather than that formula.
