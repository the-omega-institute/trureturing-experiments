# Strong grids with diagonals: independent check

This report checks the classification proposed for Problem 3 of Antony,
Lafond, Monti, and Sinaimeri, *Pairwise Compatibility Graphs of Unweighted
Caterpillars*, arXiv:2610.14681v1.

`check.py` uses only integer arithmetic. It verifies the paper's skew-grid
`[n+1,n+2]` difference labelling on several dimensions, verifies the induced
seven-vertex obstruction in the 3 by 3 king grid and its necessary-order
certificate over all 7! orders, and verifies the explicit `[0,1]` difference
labellings and `[2,3]` caterpillar witnesses for all thin grids with dimensions
at most 6 and at least two vertices. The JSON output is written by the script
next to it. Run:

```text
python3 check.py
```

The singleton `1 by 1` is omitted from the witness loop because the source's
literal degree-1 leaf convention has no finite tree with exactly one leaf.
