# HKLT Problem 4: unrestricted congruence-density bound

Source: J. Harrington, J. Klein, J. Lowrance and O. Trifonov,
*Covering systems where the prime divisors of all moduli are only 2, 3, or 5*,
arXiv:2605.18644v1, Theorem 1.9 and Section 2, Problem 4.

Run `python3 check.py`. The program uses Python's standard library, integers
and `fractions.Fraction`; its final line is `ALL_CHECKS_OK`.
The checked output is in `output.txt`.

The source fidelity checks reproduce Davenport's five-class cover of period 12
and check the source's unrestricted three-modulus assertion on every residue
assignment for three distinct moduli chosen from 2 through 8 (4,025 systems).

The counterexample has 94 distinct moduli. Its 19-class block A misses exactly
-1 modulo 5,632. The other 75 classes have residue zero and moduli m congruent
to 3 modulo 6, with 3 <= m <= 489 and 11 not dividing m. Their union is exactly
0 modulo 3 because the modulus 3 is present. Thus the whole union has period
16,896 although the least common multiple of the individual moduli is much
larger. Its uncovered residues are 11,263 and 16,895, giving density 8,447/8,448.
The checker enumerates the coprime pairs and triples independently and confirms
that the asserted upper bound is strictly smaller.

Writing b for the reciprocal sum of the added moduli, the bound is
1 - (501/5632)(b-1). The checker verifies this identity at cutoffs 489, 600 and
1,000 and an independent floor-sum certificate at 489. It checks 1,000 blocks
of the harmonic lower bound used in the unbounded family argument. The proof
of harmonic divergence is mathematical, rather than inferred from a finite run.

For the supplementary sufficient condition, a fixed-order perfect elimination
ordering of the coprimality graph makes vertices minus edges plus triangles
at least one in every nonempty induced graph. The checker verifies the pointwise
inequality for every such ordered graph on at most five vertices and every
induced subset (11,059 checks). This finite screen supports the general induction;
it does not replace it.

These are exact computational checks, not Lean kernel verification. The
preregistration stage supplies definitions and the claim only; it does not
supply a Lean proof or refutation. The checker was independently written from
the displayed arithmetic construction, without importing third-party program
or Lean source.
