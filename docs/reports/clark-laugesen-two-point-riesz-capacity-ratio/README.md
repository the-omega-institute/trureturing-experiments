# Two-point Riesz capacity ratio: exact finite checks

The target is Clark–Laugesen, *Maximizing Riesz capacity ratios: conjectures and
theorems*, arXiv:2410.14809v1, Conjecture 17. Fan restates its subrange
0 < r < s < 2 in arXiv:2609.11186v1, Conjecture 5.10.
For positive t, U_t is the maximum of the distance-power quadratic energy over
the closed probability simplex, and C_t = U_t^(1/t).

Run with Python 3, using only the standard library:

```sh
python3 docs/reports/clark-laugesen-two-point-riesz-capacity-ratio/check.py
```

`output.txt` contains the output. All arithmetic is integral or rational;
no floating point, random sampling, numerical tolerance or approximate optimizer
is used.

The fidelity checks reproduce the source's regular k-point formula
U_t = ((k-1)/k) D^t in 80 cases and its two-point mass calculation in 152 cases.
The independent check enumerates every nonzero symmetric zero-diagonal kernel
with entries in {0,1,16} on 2, 3 and 4 labels: 756 kernels in total. These include
nonmetric kernels, testing the algebraic strengthening of the proposed route.
For each kernel it computes the exact maximum at t = a/4 for
 a in {1,2,3,4,5,6,7,8,12}. The 6,804 optimized energies satisfy the diameter
lower bound and optimized-kernel comparison. All 27,216 exponent-pair comparisons
satisfy U_s^a <= 2^(b-a) U_r^b; 15,876 comparisons have 0 < r < s < 2.
Two-point kernels attain equality. Measures with zero masses are allowed.

The optimizer enumerates every nonempty support and solves the bordered
stationarity equations by exact rational Gaussian elimination. A maximum with
smallest support has a nonsingular bordered system: otherwise a zero-sum null
direction preserves the objective and permits removing a support point. The
program checks feasibility and recomputes the energy of every retained candidate.
Thus the method includes maxima on the boundary and does not assume concavity.

These results support the proposed settlement but do not prove the compact-set
statement, arbitrary positive real exponents, or absence of a prior proof.
No Lean proof of the conjecture or its negation accompanies this report.
