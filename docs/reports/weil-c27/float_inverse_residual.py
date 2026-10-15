"""Floating full-space triage, never a directed sign certificate."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
from threadpoolctl import threadpool_limits

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source-dir', type=Path, default=Path(__file__).resolve().parent)
parser.add_argument('--pilot', type=Path, default=Path(__file__).resolve().parent/'inverse-residual-pilot.json.zlib')
parser.add_argument('--out', type=Path, required=True)
args = parser.parse_args()
sys.path.insert(0, str(args.source_dir.resolve()))
from arb_far_coupling import pinned_packet, reconstruct_factors
from arb_middle_scalars import published_rows
from arb_inverse_residual_pilot import MIDDLE_SHA256, FAR_SHA256

if args.out.exists():
    raise ValueError('Use a fresh output path.')
started = time.monotonic()
scalar, rows = published_rows(args.source_dir/'arithmetic-scalars.json.zlib',
                              args.source_dir/'diagonal-floor-witness.json')
middle = pinned_packet(args.source_dir/'middle-scalars.json.zlib', MIDDLE_SHA256)
rows.update({row['n']: row for row in middle['rows']})
pilot = pinned_packet(args.pilot, '6db5093bf78fe1dcf43cbbe1b39ef2b1faf4174cf7cf33fe796c0d10ad4d231b')
if (pilot['N'], pilot['M'], pilot['middle_packet_sha256'], pilot['far_packet_sha256']) != (1280, 32768, MIDDLE_SHA256, FAR_SHA256):
    raise ValueError('Use the declared common c27 inputs.')
if not set(range(32769)) <= set(rows):
    raise ValueError('Incomplete low or middle coverage.')
far = pinned_packet(args.source_dir/'far-coupling.json.zlib', FAR_SHA256)
bits = scalar['grid_bits']
def scalar_value(row, key):
    return math.ldexp(int(row[key][0])+int(row[key][1]), -bits-1)
s = np.array([scalar_value(rows[n], 's') for n in range(32769)])
a = np.array([scalar_value(rows[n], 'a') for n in range(32769)])
N, added = 1280, 128
last = N+added
indices = np.arange(last+1)
symbols = indices*s[:last+1]
phases = np.where(indices % 2, -1., 1.)
normalization = np.ones(last+1)
normalization[0] = 1/math.sqrt(2)
denominator = indices[:, None]**2-indices[None, :]**2
actual = np.divide(symbols[None, :]-symbols[:, None], denominator,
                   out=np.zeros_like(denominator, dtype=float), where=denominator != 0)
actual *= 2/math.pi*(phases*normalization)[:, None]*(phases*normalization)[None, :]
diagonal = a[:last+1].copy()
diagonal[1:] -= s[1:last+1]/(math.pi*indices[1:])
np.fill_diagonal(actual, diagonal)
G = np.zeros_like(actual)
weights = np.ldexp(np.array([float(int(x)) for x in pilot['middle_weight_numerators']]),
                   -pilot['middle_weight_bits'])

with threadpool_limits(limits=1):
    for first in range(N+1, 32769, 1024):
        modes = np.arange(first, min(first+1024, 32769))
        den = modes[:, None]**2-indices[None, :]**2
        block = np.divide(symbols[None, :]-modes[:, None]*s[modes, None], den,
                          out=np.zeros_like(den, dtype=float), where=den != 0)
        block *= 2/math.pi*np.where(modes % 2, -1., 1.)[:, None]*(phases*normalization)[None, :]
        for row, m in enumerate(modes):
            if m <= last:
                block[row, m] = diagonal[m]
        block *= np.sqrt(weights[modes-N-1])[:, None]
        G += block.T@block
    print('middle_gram_seconds', time.monotonic()-started, flush=True)
    extension = dict(far, N=last)
    extended_scalar = dict(scalar, rows=[rows[n] for n in range(last+1)])
    alpha, D, P0, H0, _ = reconstruct_factors(extension, extended_scalar)
    P = np.array([[float(P0[i, j]) for j in range(last+1)] for i in range(P0.nrows())])
    H = np.array([[float(H0[i, j]) for j in range(H0.ncols())] for i in range(H0.nrows())])
    G += float(alpha*D)*(P.T@H@P)
    A = actual[:N+1, :N+1]
    old = A-G[:N+1, :N+1]
    cases = []
    for count in (64, 128):
        if count:
            F = actual[N+1:N+1+count, N+1:N+1+count]
            B = actual[N+1:N+1+count, :N+1]
            K = G[N+1:N+1+count, N+1:N+1+count]-F
            T = G[N+1:N+1+count, :N+1]-B
            X = np.linalg.solve(K, T)
            target = old+T.T@X
        else:
            target = old
        target = (target+target.T)/2
        values, vectors = np.linalg.eigh(target)
        minimum_vector = vectors[:, 0]
        candidate = minimum_vector.copy()
        candidate[128:] = 0
        candidate = np.rint(candidate*(1 << 40))/(1 << 40)
        candidate_trial = X@candidate
        order = np.argsort(np.abs(minimum_vector))[::-1]
        case = {'columns': count, 'smallest20': values[:20].tolist(),
                'negative_below_minus_1e_7': int(np.count_nonzero(values < -1e-7)),
                'near_zero_1e_7': int(np.count_nonzero(np.abs(values) <= 1e-7)),
                'seconds': time.monotonic()-started,
                'condition_K': float(np.linalg.cond(K)),
                'top32_vector_modes': [[int(n), float(minimum_vector[n])] for n in order[:32]],
                'candidate_dyadic_bits': 40,
                'truncated_candidate_last_mode': 127,
                'truncated_candidate_quadratic': float(candidate@target@candidate),
                'truncated_candidate_trial_integers': [str(round(float(t)*(1 << 40))) for t in candidate_trial],
                'candidate_dyadic_integers': [str(round(float(t)*(1 << 40))) for t in candidate[:128]]}
        cases.append(case)
        print(json.dumps({k:v for k,v in case.items() if k not in ('candidate_dyadic_integers', 'truncated_candidate_trial_integers')}), flush=True)
result = {'scope': 'Floating proposal only: actual finite matrix and full middle Gram plus the reused far midpoint. No far, scalar or floating error allowance, optimality verification, sign certificate, Lean or RH result.',
          'N': N, 'added': added, 'numpy': np.__version__, 'cases': cases,
          'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'input_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
              [args.pilot, args.source_dir/'arithmetic-scalars.json.zlib', args.source_dir/'diagonal-floor-witness.json', args.source_dir/'middle-scalars.json.zlib', args.source_dir/'far-coupling.json.zlib']}}
args.out.write_text(json.dumps(result, indent=2)+'\n')
