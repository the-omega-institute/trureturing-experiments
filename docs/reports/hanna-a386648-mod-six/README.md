# Hanna A386648 modulo six

Source: Paul D. Hanna, [OEIS A386648](https://oeis.org/A386648), revision 9, first COMMENTS conjecture and NAME/FORMULA (1). The 20 displayed sequence values and the Lambert W coefficient expansion used as source data retain OEIS attribution and CC BY-SA 4.0 terms. The verification code is Apache-2.0 under the repository LICENSE.

Run from the repository root:

```sh
python3 docs/reports/hanna-a386648-mod-six/check.py
```

`check.py` uses Python integers exclusively. Lists store exponential coefficients n! [x^n] F, and multiplication uses binomial convolution. Every factorial division checks exact divisibility. The implicit coefficients are generated from the original sum (A^m + LambertW(-x))^m/m!, without using the predicted congruence.

`output.txt` records a successful run. It reproduces all 20 displayed source values, verifies both source equations through degree 36, and checks all 30 conjectural residues at indices 7 through 36. It also checks the proposed tail reduction and explicit congruence through degree 36. The elementary closed expression alone is tested modulo six through index 10000; the implicit series is independently computed only through degree 36.

The final line is `ALL_CHECKS_OK`. These finite exact checks support a preregistration draft. They establish neither an all-index proof nor a Lean settlement.
