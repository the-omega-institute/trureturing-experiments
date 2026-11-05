# Hanna A357398: period four modulo ten

The target is Paul D. Hanna's first COMMENTS conjecture in [OEIS A357398](https://oeis.org/A357398), revision 13: after index zero the last digits repeat `[1,3,7,3]`. The related integral ordinary series is [A357397](https://oeis.org/A357397), revision 13.

Run with Python 3, using only its standard library:

```sh
python3 docs/reports/mod-ten-period-four/check.py
```

`output.txt` contains the result. Integer and `fractions.Fraction` arithmetic independently solve the original exponential-series equation through degree 16 and the ordinary-series equation through degree 32. The program reproduces the published A357398 terms at indices 0–16 and A357397 terms at indices 0–19, verifies the ordinary-series four-jet and its summands, compares direct substitution with the Stirling transform, and checks the proposed residues at every index 1–32.

These are finite exact checks supporting source fidelity and a proposed proof route. They do not prove the all-index conjecture and contain no Lean proof. Source sequence values are quoted for identification and fidelity under the [OEIS End-User License Agreement](https://oeis.org/LICENSE).
