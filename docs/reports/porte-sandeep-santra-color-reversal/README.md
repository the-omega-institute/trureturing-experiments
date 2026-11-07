# Porte–Sandeep–Santra color-reversal check

`check.py` implements the source local inversion literally on the labelled graph (P_4), using integer edge and sign masks. It exhaustively breadth-first searches the finite state graph from the all-positive coloring, verifies the source involution (B_{aa}=B), and checks exact shortest lengths:

- global reversal: 8, attained by `[1, 2, 3, 0, 3, 0, 1, 2]`;
- reversal of endpoints 0 and 3: 10, attained by `[0, 1, 0, 2, 0, 1, 3, 2, 3, 2]`.

The state space has eleven reachable graph masks and sixteen sign masks. Since local inversions stay inside disconnected components, projection gives additive lower bounds for (k) disjoint copies of (P_4), while concatenation gives matching upper bounds (8k) and (10k). The checker uses no floating point and ends with `ALL_CHECKS_OK`.

Run from this directory with:

```text
python3 check.py
```
