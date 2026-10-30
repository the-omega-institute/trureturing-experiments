# Semimonotone matrices: exact-order negative-entry checks

Source: Bharat Pratap Chauhan and Dipti Dubey, *On semimonotone matrices of exact order two*, arXiv:2603.00639v1, Example 3.2 and Remark 4.1 (first bullet). The question asks whether every row and column of an E_{0_k} or E_k matrix has at least k strictly negative entries, for n >= 3 and 1 <= k < n.

Run with Python 3.9 or later; the standard library suffices:

```sh
python3 docs/reports/semimonotone-negative-entries/check.py
```

`check-output.txt` contains the output. Successful execution exits 0 and ends with `ALL_CHECKS_OK`.

The checker decides the defining universally quantified matrix predicates by considering every nonempty support. On a support S, failure is equivalent to a feasible rational system x_i >= 1 and Bx <= -1 for ordinary semimonotonicity, or Bx <= 0 for strict semimonotonicity. Homogeneity justifies scaling all strictly positive quantities to at least 1. Integer Fourier--Motzkin elimination decides feasibility, and exact Fraction back substitution checks each feasible witness against every original inequality. No floating point is used. All 81 order-two matrices over {-1,0,1} are also checked against the independent determinant criterion of Tsatsomeros--Wendler (2019), Propositions 3.6--3.7; zero-boundary controls distinguish the two predicates.

Fidelity: Example 3.2 has [4,6,0,0] semimonotone principal blocks of sizes 1 through 4, respectively, and exact order 2. Its displayed inverse is reproduced and its determinant is -5/4, all rationally.

Finite screen: all 19,683 matrices of order 3 over {-1,0,1}, and all 4,096 matrices of order 4 with diagonal 1 and off-diagonal entries in {-1,0}. There are 366 qualifying (matrix, predicate, exact-order) instances. Each satisfies the question's row and column count. Each also has a strongly connected negative off-diagonal digraph after every deletion of fewer than k vertices. Exact-order membership is decided from the matrix predicates, independently of these sign and graph tests.

These are exhaustive checks of the stated finite alphabets, not a proof over all real matrices or all dimensions. No Lean proof is included.
