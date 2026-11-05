# Beshkov ReLU sphere-width helper check

`check.py` uses only Python `fractions.Fraction` and integer arithmetic. For dimensions 1 through 4 it checks the proposed network from Conjecture 2.16 on an exact rational input grid and verifies the source's simplex-boundary image identity on exact rational boundary points at denominators 2 and 3. The inverse is the source's formula with `t = (1 - sum(y))/(2d - 1)`, including zero-coordinate faces. It also checks the dimensionally consistent interpretation of the source's sums: the last hidden preactivation sums the `d` actual inputs, equivalently appending a fixed coordinate `x_(d+1)=0`.

These are helper checks for the planned all-dimensional construction and do not prove the topological homeomorphism or the lower-bound obstruction. No floating-point values are used as evidence. A successful run ends with `ALL_CHECKS_OK`.
