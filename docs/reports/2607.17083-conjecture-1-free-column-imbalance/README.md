# Free-column imbalance in constrained binary covering arrays

Exact preregistration checks for Conjecture 1, “Balanced Columns for Unconstrained Parameters”, in Serra et al., *Optimal Combinatorial Testing with Constraints: The Balancing Act*, arXiv:2607.17083v1, Section 4.4.

Run with Python 3 and `z3-solver==4.15.4`:

```sh
python3 check.py
```

The fidelity check copies the source’s matrix M^(1), verifies all 180 pairwise interactions and balanced columns, and reproduces CAN(2,10,2)=6 using the source’s binomial formula.

For the proposed family, forbid (x_i,y_j)=(1,1) for i unequal to j and leave z free. The program checks the displayed b+2-row cover for 3 ≤ b ≤ 9. It enumerates all complete assignments at b=1,...,5 to verify reducedness (every symbol occurs and exactly the listed pairs are forbidden).

The optimum checks encode the covering problem as nonnegative integer multiplicities of every valid complete test. Each feasible interaction must be covered at least once. Z3 solves these integer linear constraints in QF_LIA; all coefficients, bounds, and outputs are exact integers. Duplicate tests are allowed. The program excludes every smaller row count for b=3,4,5, exhausts all optimal row multisets, and checks every possible frequency of z. There are exactly two optimal row multisets, exchanged by complementing z. At b=4 the optimum is 6 and z has counts 2 and 4; at b=5 the optimum is 7 and counts are 2 and 5. For b=4 a balanced z requires exactly 7 rows.

Boundary checks give optima 4, 5 and 5 for b=1, 2 and 3, respectively, with balanced optima. The displayed b+2-row construction itself fails at b=1,2. For b=3 its counts 2 and 3 are balanced.

`output.txt` contains the program’s output, ending in `ALL_CHECKS_OK`. These are exact external computations; the general family argument and the negation of the source conjecture remain to be checked in Lean after preregistration. No Lean proof is included here.
