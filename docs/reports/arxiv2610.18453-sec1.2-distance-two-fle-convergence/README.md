# Distance-two FLE convergence check

This report checks a structural identity used in Armfield and Mailler, *The forward-loop-erased ants walk*, arXiv:2610.18453v1, Section 1.2 and Theorem 1.1. For every enumerated finite simple connected graph on three, four, or five vertices with nest–food distance two, and two positive integer weight assignments, `check.py` solves the absorbing random-walk equations over exact rational numbers. It verifies, for every common neighbour `i`,

`p_{N i}/W_{N i} + p_{i F}/W_{i F} = 1/C_{N,F}`.

The script also performs integer-draw simulations of the forward-loop-erased update on a 4-cycle with opposite nest and food and on `K_{2,3}` for 5,000 ants. The normalized terminal-edge weights are printed at 100, 1,000, and 5,000 ants as a sanity reading. The simulation is not used as proof evidence.

Run:

```sh
python3 check.py | tee check.out
```

The final line is `ALL_CHECKS_OK` and the process exits 0.
