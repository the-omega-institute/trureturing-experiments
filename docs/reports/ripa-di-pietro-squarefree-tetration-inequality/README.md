# Composite squarefree tetration-speed inequality

The target is the last branch of equation (3.4) in Appendix B, Conjecture 3.1 of Marco Ripà and Gabriele Di Pietro, *A compact notation for peculiar properties characterizing integer tetration*, NNTDM 32(3) (2026), 477–496, DOI [10.7546/nntdm.2026.32.3.477-496](https://doi.org/10.7546/nntdm.2026.32.3.477-496). The final author source is [Zenodo v12, record 23404787](https://zenodo.org/records/23404787), licensed CC BY 4.0.

For every composite squarefree radix r and a > 1 with r not dividing a, the proposed inequality compares the eventual radix-r speed of a with the minimum eventual speed of its distinct prime divisors. Noncoprime pairs are included.

Run with Python 3, with no third-party dependencies:

    python3 check.py

`check.out` contains the complete output. The program uses exact integers exclusively. Modular tetration is evaluated separately on prime-power factors: Euler reduction for unit bases, exact exponent thresholds for nonunits, and CRT. Direct integer towers cross-check the evaluator. It reproduces all seven journal residues for base 27057 modulo 10^25.

For every composite squarefree radix at most 120 and every base from 2 through 150 not divisible by that radix, it checks the proposed local-slope formula at heights 8 through 12 and checks the minimum-factor inequality. It also checks the source’s established prime-radix formulas. Orders truncated at precision 128 are treated as lower bounds; each finite minimum is required to be strictly below that precision, so every reported radix order is exact. This includes noncoprime bases, and avoids treating exponent reduction as valid on nonunits.

Finite agreement does not prove eventual constancy, identify the eventual value for every input, or prove the unbounded conjecture. The intended settlement uses a general argument through stabilization, LTE, nonunit growth, eventual minima of affine functions, and multiplicative congruences; enumeration is supporting evidence only.
