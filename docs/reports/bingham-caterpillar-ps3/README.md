# Caterpillar principal specialization family

Source: Aram Bingham, Kyle Celano, Patricia Commins, Rosa Orellana, Jianping Pan,
and Mary Claire Simone, *Principal specializations of the t-chromatic symmetric
function*, arXiv:2610.31785v1, Remark 3.11, first question.

Run with Python 3 and SymPy 1.14.0:

```sh
python3 docs/reports/bingham-caterpillar-ps3/check.py
```

The stored output is `output.txt`. All arithmetic is exact. The program constructs
the source's explicit caterpillar vertex and edge sets and independently enumerates
proper colorings for the base pair (5,4,3)/(4,5,3) and both members of the family
(5,j,j,4,3)/(4,j,j,5,3), for j=1,2,3. These counts agree with the separate leaf-summed
transfer calculation. It reproduces the source's reported j=1..7 behavior and
checks j=8..30. Breadth-first bipartitions independently check the two-color formulas.

Finally, the program expands the proposed arbitrary-diagonal matrix identity in
Z[q,x,y,z] and obtains the zero polynomial. This is exact algebraic evidence for a
general proof route. Finite family checks alone do not prove the universal claim;
the required Lean theorem and combinatorial transfer-factorization proof remain
outside this experiment. The report concerns the ordinary chromatic symmetric
function X, with no color-count variable t, and does not address the separate
classification question or either four-color conjecture.

Expected last line: `ALL_CHECKS_OK`.
