# Del Valle–Praeger digraph-code lift

`check.py` is an exact GF(2) finite consistency check for the proposed
three-layer/two-fibre lift of a graph code.  It uses the complete bipartite
graph K₃,₃ and the even-parity [3,2,2] local code, reproducing the graph-code
parameters [9,4,4] and checking the lifted 216-coordinate code, its rank,
minimum distance, properness, strong connectivity, unequal input dimensions,
and nine finite symmetry generators.  All arithmetic is integer or exact
bit-vector arithmetic; no floating-point values are used.

Run:

```sh
python3 check.py
```

The recorded output is in `check.out`; success ends with `ALL_CHECKS_OK`.
