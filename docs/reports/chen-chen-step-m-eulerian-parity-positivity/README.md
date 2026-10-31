# Chen–Chen even-row positivity checks

The targets are the three conjectures in [OEIS A392664](https://oeis.org/A392664/internal), [A392665](https://oeis.org/A392665/internal), and [A392666](https://oeis.org/A392666/internal), for fixed parameters m = −3, −4, and −5. Their endpoints equal 1; their interior recurrence is

T(n,k) = (m(n−k)+1) T(n−1,k−1) + (m(k−1)+1) T(n−1,k).

Run with Python 3 and its standard library:

    python3 docs/reports/chen-chen-step-m-eulerian-parity-positivity/check.py

The stored output is `check.out`. All arithmetic uses exact Python integers. The script reproduces all 165 entries in the displayed rows 1–10 of the three sources and checks them against an independent memoized evaluation. It reproduces the sources' verification through row 400: 120,600 even-row entries are strictly positive. It also checks reflection, the division-free row-sum relation, even-row alternating sums, and the proposed two-step second-column identity and lower bounds. A separate grid checks parameters m = −3 through −32 and rows 1–80.

These finite checks support an induction route. They do not prove positivity for all even rows or for all parameters. The proposed unbounded proof handles the second and penultimate columns separately, then propagates signs across pairs of rows. No Lean proof is included.

The displayed row data and recurrence are mathematical facts transcribed from the three OEIS entries. The related summary is [Chen and Chen, Zenodo v2, DOI 10.5281/zenodo.18452218](https://zenodo.org/records/18452218); the registered target is the conjunction of the three OEIS conjectures.
