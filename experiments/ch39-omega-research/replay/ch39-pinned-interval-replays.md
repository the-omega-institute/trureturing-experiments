# Replaying two interval certificates

Both replays use Python with `python-flint==0.9.0`. Run with assertions enabled. Supply the replay directory from the indicated input commit. This fixes the shell geometry, operators, targets, and other numerical data consumed by each script.

## Finite-height exclusion for the ray readout

The input commit is `d0c5ab50ebaf763fea12e2629cb6760205448f08`. Run both grids with the same pinned input.

```sh
python ch39-ray-blind-height-scan.py "$ray_replay" --grid-bits 180
python ch39-ray-blind-height-scan.py "$ray_replay" --grid-bits 200
```

Both runs cover the closed height interval from $14$ to $14\cdot 2^{80}$. On each dyadic subinterval, the script encloses the candidate blind phases for all five charges and all three cyclic characters. An exact integer residue count then excludes hits on the clock trajectory. The first positive zeta zero is enclosed above $14$. [Derivation and original script](https://github.com/the-omega-institute/trureturing-books/issues/1687#issuecomment-7233848461)

## Same-sign repeated gate with one free phase

The input commit is `9e2b8abfe862acb15667d5a0154731a67e491a95`.

```sh
python ch39-same-sign-phase-circle.py "$phase_replay"
```

The replay fixes eleven angles at the decimal center of the sign-flipped certificate. It varies the common phase of the repeated 52-chord gate around the full circle and checks the six residuals for shells seven through nine. In 316 interval evaluations the terminal exclusion counts are `[123, 3, 20, 0, 28, 0]`, with zero unresolved arcs. The feasibility of a same-sign solution with all gate angles free remains open. [Derivation and original script](https://github.com/the-omega-institute/trureturing-books/issues/1627#issuecomment-7234066551)
