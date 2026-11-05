# Goetze–Zang Hypothesis 19: chain/crown separation

`check.py` is an independent Python 3 check for arXiv:2609.09873v1, Hypothesis 19. It uses integer edge masks, finite set enumeration, and `fractions.Fraction`; no floating-point value is used as evidence.

The script verifies the crown family `H_n` has the two ordinary chain layers `i < j` and `i > j`, rejects a one-layer cover by exhaustive induced-`2K2` enumeration, and exhaustively enumerates weak-induced chain layers for `n = 3, 4, 5`. The exact readings are `ic_u(H_3) = ic_u(H_4) = 2` and `ic_u(H_5) = 3`; restriction to the induced `H_5` subhost gives the lower bound for `n = 6, 7, 8`, while binary labels give the matching upper bound 3. It also checks the component intersection lemma on all finite support pairs through `n = 8` and prints the family parameters `n = 8^k + 1`.

Run:

```text
python3 check.py
```

The stored output ends with `ALL_CHECKS_OK`.
