# Exact checks for stabilizer-sector non-T₀ candidate

Source: Cathy Cheng, Nigel Higson, Danilo Polo Ojito and Emil Prodan,
*Ground-state spaces of generalized stabilizer codes*, arXiv:2610.33221v1,
Conjecture 5.6, equations (4.7) and (5.1), and Proposition 6.23.
Source text: https://arxiv.org/abs/2610.33221v1 .

Run with Python 3, using only the standard library:

```sh
python3 docs/reports/chpp-2610-stabilizer-sector-t0/check.py
```

The stored result is `output.txt`. The final line is `ALL_CHECKS_OK`.
All calculations use integers, finite sets and binary parity. A spin bit 0
means +1 and bit 1 means −1.

The fidelity check reproduces the biased one-dimensional Ising boundary
criterion of Proposition 6.23 on all 512 configurations with a seven-site
variable window and arbitrary constant tails, testing all 512 nine-site flips
for each. Exactly ten boundary configurations occur in this window, with four
asymptotic sector labels. It also checks a movable zero-energy domain wall.

The candidate is the diagonal qubit model on the connected comb in ℤ²:
one origin pin, every horizontal nearest-neighbor edge and the vertical axis
edges. The checks verify 81 finite path-generation identities, the incidence
bound of five generators per site, and the symplectic/Lagrangian calculation
on three sites. They then check 110,592 finite spin flips of row-wall
configurations, comparing equation (4.7) to an independently counted local
Hamiltonian difference and to the opposite-tail row-wall lower bound.
There are 59,049 exact finite wall-shift pairs and 4,096 observation-cylinder
checks in both directions between the constant wall positions 0 and 1.

Each configuration is specified on the infinite lattice; every energy change
is computed using the finitely affected generators. No divergent total energy
is subtracted, and no truncated Hamiltonian is substituted for the model.
The finite checks support the proposed construction. They do not prove the
infinite ground-configuration assertion, the Pontryagin-dual identification,
or failure of T₀ for the infinite quotient. Those require the general argument
and a subsequent Lean probe; no Lean settlement is supplied by this report.

The checker is original code and introduces no third-party program dependency.
