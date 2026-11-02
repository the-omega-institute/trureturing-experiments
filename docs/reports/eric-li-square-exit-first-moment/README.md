# Exact checks for Eric Li's square-exit first moment

Run `python3 docs/reports/eric-li-square-exit-first-moment/check.py` from the repository root. The program uses only Python's standard library and exact integers and fractions. Its recorded output is `check.out`; success ends with `ALL_CHECKS_OK`.

Source: Eric Li, *Square-Annular Dynamics and Coalescence Frontiers for n + τ(n)*, arXiv:2606.17926v1, section 5's exit-set definition and exact table, Proposition 7.5, and Remark 8.25. The source manuscript is dated 13 June 2026 and arXiv records submission on 16 June 2026. The source text is not redistributed here.

The fidelity check reproduces every row of the source's table for levels 2 through 12, including active deficits, distinct exit offsets, and transferred exits. It also checks Proposition 7.5 through level 256. Trial factorization and an independently implemented smallest-prime-factor sieve agree through 4096.

The proposed refutation uses square deficits and finite periodic concentration. The program enumerates all 148225 residues modulo 5²·7²·11², verifying the valuation-one probabilities, independence across different primes for each fixed shift, and the exact first and second moments. It checks the square-deficit divisibility and collision implications for every eligible level at most 100000 and every J from 1 through 4. The valuation-one sufficient premise occurs 99215 times for J=1 and 1419 times for J=2. It occurs zero times for J=3 and J=4 in this window; those two tests of that sufficient premise are vacuous. The weaker common-power divisibility premise has 40355 and 3039 successful levels respectively, and the resulting active offsets are distinct in every such case.

These finite checks support the stated formulas and implications. They establish no asymptotic density, rate of growth, or refutation of a big-O bound. The planned settlement requires a general argument and a later Lean proof.
