# OEIS A388712 and A387118: partitions without choosable initial intervals

Gus Wiseman conjectures that partitions without choosable initial intervals {1..y₁}, …, {1..y_k} are equinumerous with partitions whose weakly increasing order is not superdiagonal (z_i ≥ i for all i). A388712 states this by number of parts, and A387118 states it in total.

`check.py` (Python standard library) runs over all partitions of n ≤ N. It decides choosability by searching for distinct representatives and checks that this coincides with superdiagonality pointwise. It also compares the counts with the published A388712 triangle and A387118 data.

Run it with `python3 check.py 22`. Exit 0 and `RESULT PASS` mean everything agrees. `run.log` is the output for N = 22.
