# Mugnolo Conjecture 4.7: uniform Dirichlet mean-distance gap

`check.py` performs the preregistration check with `fractions.Fraction` and SymPy exact algebra. It verifies the variance/truncation combination, the one-sided variance estimate on exact two-point distributions, and the exact ratio `λ₁ρ² = π²/16` for a path, a star, a cycle, and a figure-eight. The displayed decimal values are only high-precision readability output; the inequality uses the rational bound `π < 22/7`.

Run:

```sh
python3 check.py
```

The final line must be `ALL_CHECKS_OK`.
