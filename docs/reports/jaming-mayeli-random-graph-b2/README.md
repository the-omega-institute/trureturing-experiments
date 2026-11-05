# Jaming–Mayeli random-graph B₂ check

`check.py` uses only integer matrices and exact SymPy polynomial arithmetic. It verifies the combinatorial Laplacian definition and the spectrum/B₂ property for the three-vertex path, exhaustively checks every connected labelled graph on 2, 3, and 4 vertices, and checks four fixed connected five-vertex samples. The pair-sum polynomial is obtained from an exact resultant; squarefreeness is equivalent to distinct unordered pair sums, including repeated indices.

This is a bounded sanity check. It does not prove or refute the asymptotic Erdős–Rényi probability statement, and it uses no floating-point computation as evidence.
