# Kurkov A395726: exact checks of the three Stirling formulas

Source: Mikhail Kurkov, [OEIS A395726](https://oeis.org/A395726/internal), approved revision 21. Dependent sequences: [A002024](https://oeis.org/A002024) and [A002260](https://oeis.org/A002260). The embedded DATA and displayed rows are attributed OEIS source data; see the [OEIS license agreements](https://oeis.org/wiki/The_OEIS_End-User_License_Agreement). The checking code is independently written.

Run from the repository root with Python 3; only the standard library is required:

```sh
python3 docs/reports/kurkov-a395726-stirling-formulas/check.py
```

The expected output is in `check.out`, ending with `ALL_CHECKS_OK`. All arithmetic uses integers and `fractions.Fraction`, including the exponent −1 at n=0 in the first triangular block.

The defining recurrence reproduces all 78 published DATA terms and all 90 displayed array entries. An independent translation of the published PARI column program agrees at 6,321 pairs. Conjectures 1 and 2 are checked through row 20 and column 300, including the respective boundary cases; Conjecture 3 is checked for rows 1–20 and columns 0–300. The triangular-block Stirling identification is checked on the same positive-column rectangle.

The unnumbered source formula preceding the conjectures has an indexing discrepancy: n=1, k=2 gives 1 by NAME and 3 by its printed Stirling expression. It is outside the three numbered targets and is not a premise of any check.

Finite agreement is evidence of transcription fidelity and support for a proposed general argument. It does not prove any unbounded identity. No Lean proof is included.
