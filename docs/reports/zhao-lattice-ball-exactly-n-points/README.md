# Zhao lattice ball count check

This report checks the normalization and small planar instances of Jianqiang Zhao, *Largest Circle Enclosing Exactly n Interior Lattice Points. II*, arXiv:2609.16081v1, Conjecture 25(a). The script uses `fractions.Fraction` for all geometric distances and SymPy's exact `pi` and `gamma` expressions for the volume identity. No floating-point value is used as evidence.

For each `n = 1, ..., 12`, it searches rational centers in `[0,1]^2`, counts all lattice points in the radius-squared `n/pi` disk, and certifies every comparison with the rational bounds `333/106 < pi < 355/113`. Since `pi > 3` and `n <= 12`, enumerating `[-4,4]^2` is exhaustive. The final line is `ALL_CHECKS_OK`.
