"""Floating coefficient-uniform triage using finite displacement identities.

Two finite actual-H inverse actions reconstruct the trial on every retained
column. Finite weighted-Gram generators avoid a new dense middle Gram.
The infinite far part reuses the archived moments. No floating output is a
sign, inertia, infinite-operator, or RH certificate.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from threadpoolctl import threadpool_limits

from arb_far_coupling import reconstruct_factors
from arb_inverse_residual_obstruction import HERE, M, N, inputs


def run(args):
    started = time.monotonic()
    scalar, rows, pilot, far = inputs()
    last = args.last_mode
    if not N < last < M:
        raise ValueError('Finite selected modes must lie in the paid middle.')
    bits = scalar['grid_bits']

    def midpoint(row, key):
        return math.ldexp(int(row[key][0])+int(row[key][1]), -bits-1)

    symbol = np.array([n*midpoint(rows[n], 's') for n in range(M+1)])
    diagonal = np.array([midpoint(rows[n], 'a') for n in range(last+1)])
    diagonal[1:] -= symbol[1:last+1]/(math.pi*np.arange(1, last+1)**2)
    modes = np.arange(M+1)
    normalization = np.ones(M+1)
    normalization[0] = 1/math.sqrt(2)
    phase = np.where(modes % 2, -1., 1.)*normalization
    generators = np.column_stack((phase, phase*symbol))
    kappa = 2/math.pi
    retained = np.arange(N+1)
    selected = np.arange(N+1, last+1)
    weights = np.ldexp(np.array([float(int(x)) for x in pilot['middle_weight_numerators']]),
                       -pilot['middle_weight_bits'])

    def cross(ms, ns, left, right):
        denominator = ms[:, None]**2-ns[None, :]**2
        numerator = left[:, 0, None]*right[None, :, 1]-left[:, 1, None]*right[None, :, 0]
        return np.divide(kappa*numerator, denominator,
                         out=np.zeros_like(numerator), where=denominator != 0)

    def actual(ms, ns):
        result = cross(ms, ns, generators[ms], generators[ns])
        positions = {int(n):j for j,n in enumerate(ns)}
        for i,m in enumerate(ms):
            if int(m) in positions:
                result[i, positions[int(m)]] = diagonal[m]
        return result

    with threadpool_limits(limits=1):
        F = actual(selected, selected)
        B = actual(selected, retained)
        U = np.linalg.solve(F, generators[selected])
        GP = generators[retained]-B.T@U
        Y = cross(selected, retained, U, GP)
        schur = cross(retained, retained, GP, GP)
        np.fill_diagonal(schur, diagonal[:N+1]-np.sum(B*Y, axis=0))
        solve_residual = B-F@Y
        inverse_residual_max = float(np.max(np.abs(solve_residual)))
        inverse_generator_residual_max = float(np.max(np.abs(F@U-generators[selected])))
        print('two-RHS trial reconstructed; seconds', time.monotonic()-started, flush=True)

        T = np.zeros((N+1, 2))
        gram_diagonal = np.zeros(N+1)
        checks = np.unique(np.array([0,1,5,32,64,128,512,1024,N-1,N]))
        direct_check_gram = np.zeros((len(checks), len(checks)))
        for first in range(last+1, M+1, 512):
            ms = np.arange(first, min(first+512, M+1))
            GO = generators[ms]-actual(ms, selected)@U
            residual = cross(ms, retained, GO, GP)
            w = weights[ms-N-1]
            T += residual.T@(w[:, None]*GO)
            gram_diagonal += np.sum(w[:, None]*residual**2, axis=0)
            check_rows = residual[:, checks]
            direct_check_gram += check_rows.T@(w[:, None]*check_rows)

        # These are finite sums only. No infinite generator sum is taken.
        gram = cross(retained, retained, T, GP)+cross(retained, retained, GP, T)
        np.fill_diagonal(gram, gram_diagonal)
        finite_target = schur-gram
        middle_check_max = float(np.max(np.abs(gram[np.ix_(checks,checks)]-direct_check_gram)))
        print('finite middle generators accumulated; seconds', time.monotonic()-started, flush=True)

        alpha, scale, Pq, Hq, _ = reconstruct_factors(
            dict(far, N=last), dict(scalar, rows=[rows[n] for n in range(last+1)]))
        P = np.array([[float(Pq[i,j]) for j in range(last+1)] for i in range(Pq.nrows())])
        H = np.array([[float(Hq[i,j]) for j in range(Hq.ncols())] for i in range(Hq.nrows())])
        feature = P[:, :N+1]-P[:, N+1:]@Y
        far_gram = float(alpha*scale)*(feature.T@H@feature)
        target = finite_target-far_gram
        symmetry_error = float(np.max(np.abs(target-target.T)))
        target = (target+target.T)/2
        values, vectors = np.linalg.eigh(target)
        v = vectors[:, 0]
        dyadic_bits = 64
        dyadic = [str(round(float(t)*(1 << dyadic_bits))) for t in v]
        test_vector = np.ldexp(np.array([float(int(t)) for t in dyadic]), -dyadic_bits)
        finite_inverse_trial = Y@test_vector
        result = {
            'schema':'ACTUAL_C27_STRUCTURED_JOINT_FLOAT_V1',
            'N':N, 'M':M, 'last_mode':last,
            'basis':scalar['basis'],
            'algorithm':'Two RHS for F^-1 G_E; finite displacement Schur/residual/weighted-Gram reconstruction; archived far feature transform.',
            'finite_generator_rhs_count':2,
            'smallest32':values[:32].tolist(),
            'negative_below_minus_1e_7':int(np.count_nonzero(values < -1e-7)),
            'near_zero_1e_7':int(np.count_nonzero(np.abs(values) <= 1e-7)),
            'minimum_candidate_dyadic_bits':dyadic_bits,
            'minimum_candidate_dyadic_integers':dyadic,
            'minimum_candidate_midpoint_quadratic':float(test_vector@target@test_vector),
            'minimum_candidate_finite_Q':float(test_vector@schur@test_vector),
            'minimum_candidate_middle_debit':float(test_vector@gram@test_vector),
            'minimum_candidate_far_debit':float(test_vector@far_gram@test_vector),
            'minimum_candidate_trial_dyadic_bits':dyadic_bits,
            'minimum_candidate_trial_dyadic_integers':[str(round(float(t)*(1 << dyadic_bits))) for t in finite_inverse_trial],
            'diagnostics':{
                'finite_inverse_generator_residual_max':inverse_generator_residual_max,
                'reconstructed_all_rhs_residual_max':inverse_residual_max,
                'middle_direct_subset_modes':checks.tolist(),
                'middle_direct_subset_gram_max_difference':middle_check_max,
                'target_symmetry_error':symmetry_error,
            },
            'elapsed_seconds':time.monotonic()-started,
            'numpy':np.__version__,
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'dependency_sha256':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
                ['arb_arithmetic.py','arb_far_coupling.py','arb_tail_moments.py','arb_middle_scalars.py',
                 'arb_inverse_residual_pilot.py','arb_inverse_residual_obstruction.py','diagonal_floor_witness.py']},
            'input_sha256':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
                ['arithmetic-scalars.json.zlib','middle-scalars.json.zlib','diagonal-floor-witness.json',
                 'far-coupling.json.zlib','inverse-residual-pilot.json.zlib']},
            'scope':'Floating proposal only on the full retained coefficient space. No scalar, moment or old Gram production rerun. '
                    'No directed arithmetic, inverse, finite-rounding or far allowance is included. '
                    'Finite selected residual is zero in the exact two-RHS identities, not necessarily in floating arithmetic. '
                    'Tiny eigenvalues do not certify a sign or a subspace. No actual Schur positivity, Lean verification, cofinal positivity or RH result.',
        }
        args.out.write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({k:result[k] for k in ['smallest32','negative_below_minus_1e_7','near_zero_1e_7','diagnostics','elapsed_seconds']},indent=2),flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--last-mode', type=int, default=2048)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('Use a fresh output path.')
    run(args)


if __name__ == '__main__':
    main()
