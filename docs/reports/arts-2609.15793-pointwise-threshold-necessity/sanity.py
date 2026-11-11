#!/usr/bin/env python3
"""Floating-point joint-graph simulation; diagnostic only, never proof evidence."""
import numpy as np

rng = np.random.default_rng(260915793)
L, n, trials = 64, 4096, 5
rho = np.log(n) / n
S = np.array([[L, 2*L], [2*L, 4*L+1]])
assert rho * S.max() < 1
print('DIAGNOSTIC_ONLY: floating-point joint off-diagonal Poisson simulations')
print(f'seed=260915793 L={L} n={n} trials={trials}; rho*S_max<1')
for trial in range(trials):
    z = rng.integers(0, 2, size=n)
    theta = np.ones(n)
    for a in (0, 1):
        ids = np.flatnonzero(z == a)
        pairs = ids[:2*(len(ids)//2)]
        theta[pairs[::2]], theta[pairs[1::2]] = .99, 1.01
    degree = np.zeros(n, dtype=np.int64)
    # One common graph: each upper-triangular count contributes to both degrees.
    for i in range(n-1):
        rates = rho * theta[i] * theta[i+1:] * S[z[i], z[i+1:]]
        edges = rng.poisson(rates)
        degree[i] += edges.sum()
        degree[i+1:] += edges
    estimate = (degree >= (35*L/16)*np.log(n)).astype(np.int64)
    errors = np.count_nonzero(estimate != z)
    print(f'trial={trial+1} block_sizes=({np.count_nonzero(z==0)},{np.count_nonzero(z==1)}) errors={errors}')
print('DIAGNOSTIC_DONE')
