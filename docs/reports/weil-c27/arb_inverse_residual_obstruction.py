"""Directed tests of fixed-column c27 inverse-residual certificates.

Reuse the complete scalar packets, archived weights and far moments. The
input direction is the existing full1281-coordinate floating proposal,
rounded to exact dyadics. No new eigenvector or scalar table is generated.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time

import flint
from flint import arb, arb_mat, ctx, fmpq

from arb_arithmetic import Arithmetic, grid_interval
from arb_far_coupling import far_bounds, pinned_packet, reconstruct_factors
from arb_inverse_residual_pilot import MIDDLE_SHA256, FAR_SHA256, trial_integer
from arb_middle_scalars import published_rows
from diagonal_floor_witness import interval_ball


HERE = Path(__file__).resolve().parent
PILOT_SHA = '6db5093bf78fe1dcf43cbbe1b39ef2b1faf4174cf7cf33fe796c0d10ad4d231b'
N, M = 1280, 32768


def inputs():
    scalar, rows = published_rows(HERE/'arithmetic-scalars.json.zlib',
                                  HERE/'diagonal-floor-witness.json')
    middle = pinned_packet(HERE/'middle-scalars.json.zlib', MIDDLE_SHA256)
    if set(rows) & {r['n'] for r in middle['rows']}:
        raise ValueError('Overlapping scalar coverage.')
    rows.update({r['n']: r for r in middle['rows']})
    if set(rows) != set(range(M+2)):
        raise ValueError('Require the complete published scalar coverage.')
    pilot = pinned_packet(HERE/'inverse-residual-pilot.json.zlib', PILOT_SHA)
    far = pinned_packet(HERE/'far-coupling.json.zlib', FAR_SHA256)
    if (far['N'], far['M'], far['jet_order'], far['moment_diagonals'],
            pilot['N'], pilot['M']) != (N, M, 56, 56, N, M):
        raise ValueError('Wrong common-form parameters.')
    return scalar, rows, pilot, far


def direction(path):
    candidate = json.loads(path.read_text())
    if (candidate['N'], len(candidate['candidate_dyadic_integers'])) != (N, N+1) or candidate['columns'] not in (64,128):
        raise ValueError('Require the untruncated candidate in the original mass.')
    return candidate


def propose(args):
    import numpy as np
    from threadpoolctl import threadpool_limits
    scalar, rows, pilot, far = inputs()
    candidate = direction(args.candidate)
    count = candidate['columns']
    selected = tuple(range(N+1,N+count+1))
    bits = scalar['grid_bits']
    def mid(row, name):
        return math.ldexp(int(row[name][0])+int(row[name][1]), -bits-1)
    s = np.array([mid(rows[n], 's') for n in range(M+1)])
    a = np.array([mid(rows[n], 'a') for n in range(M+1)])
    last = selected[-1]
    indices = np.arange(last+1)
    phase = np.where(indices % 2, -1., 1.)
    norm = np.ones(last+1)
    norm[0] = 1/math.sqrt(2)
    symbol = indices*s[:last+1]
    diagonal = a[:last+1].copy()
    diagonal[1:] -= s[1:last+1]/(math.pi*indices[1:])
    def block(modes):
        den = modes[:, None]**2-indices[None, :]**2
        b = np.divide(symbol[None, :]-modes[:, None]*s[modes, None], den,
                      out=np.zeros_like(den, dtype=float), where=den != 0)
        row_norm = np.where(modes == 0, 1/math.sqrt(2), 1.)
        b *= (2/math.pi*np.where(modes % 2, -1., 1.)*row_norm)[:, None]*(phase*norm)[None, :]
        for i, m in enumerate(modes):
            if m <= last:
                b[i, m] = diagonal[m]
        return b
    v = np.ldexp(np.array([float(int(x)) for x in candidate['candidate_dyadic_integers']]),
                 -candidate['candidate_dyadic_bits'])
    weights = np.ldexp(np.array([float(int(x)) for x in pilot['middle_weight_numerators']]),
                      -pilot['middle_weight_bits'])
    joint = np.zeros((count+1, count+1))
    with threadpool_limits(limits=1):
        selected_rows = block(np.array(selected))
        F, Bv = selected_rows[:, N+1:], selected_rows[:, :N+1]@v
        for first in range(N+1, M+1, 1024):
            modes = np.arange(first, min(first+1024, M+1))
            b = block(modes)
            cols = np.column_stack((b[:, :N+1]@v, b[:, N+1:]))
            joint += cols.T@(weights[modes-N-1, None]*cols)
        alpha, D, P0, H0, _ = reconstruct_factors(
            dict(far, N=last), dict(scalar, rows=[rows[n] for n in range(last+1)]))
        P = np.array([[float(P0[i,j]) for j in range(last+1)] for i in range(P0.nrows())])
        H = np.array([[float(H0[i,j]) for j in range(H0.ncols())] for i in range(H0.nrows())])
        features = np.column_stack((P[:, :N+1]@v, P[:, N+1:]))
        joint += float(alpha*D)*(features.T@H@features)
        K, t = joint[1:, 1:]-F, joint[1:, 0]-Bv
        x = np.linalg.solve(K, t)
        partial = selected_rows[:, N+1:]
        partial_K = partial.T@(weights[:count, None]*partial)-F
        print('partial_K_min', np.linalg.eigvalsh(partial_K)[0], flush=True)
        print('full_K_min', np.linalg.eigvalsh(K)[0], flush=True)
        retained_rows = block(np.arange(N+1))
        q = float(v@retained_rows[:, :N+1]@v)
        margin = q-joint[0,0]+2*x@t-x@K@x
        print('fixed_direction_proposed_margin', margin, flush=True)
    candidate['trial_dyadic_bits'] = 64
    candidate['trial_dyadic_integers'] = [str(round(float(t)*(1 << 64))) for t in x]
    candidate['proposal_scope'] = 'Midpoint solve for this fixed direction only; no interval conclusion.'
    args.out.write_text(json.dumps(candidate, indent=2)+'\n')


def certify(args):
    if args.output_bits < 64 or args.precision < max(640,args.output_bits)+96:
        raise ValueError('The reused640-bit far supplier requires at least96 guard bits.')
    started = time.monotonic()
    scalar, rows, pilot, far = inputs()
    candidate = direction(args.candidate)
    count = candidate['columns']
    selected = tuple(range(N+1,N+count+1))
    if len(candidate['trial_dyadic_integers']) != count:
        raise ValueError('Require the complete selected trial.')
    vq = [fmpq(int(i), 1 << candidate['candidate_dyadic_bits'])
          for i in candidate['candidate_dyadic_integers']]
    xq = [fmpq(int(i), 1 << candidate['trial_dyadic_bits'])
          for i in candidate['trial_dyadic_integers']]
    uq = vq+[-x for x in xq]
    u = [arb(x) for x in uq]
    ucolumn = arb_mat([[x] for x in u])
    norm_u_sq = sum((x*x for x in uq), fmpq(0))
    norm_v_sq = sum((x*x for x in vq), fmpq(0))
    if not norm_v_sq > 0:
        raise ValueError('Nonzero retained witness required.')
    data = {n: (interval_ball(row['s'], scalar['grid_bits']),
                interval_ball(row['a'], scalar['grid_bits']))
            for n, row in rows.items()}
    arithmetic = Arithmetic(27, scalar['grid_bits'])
    weights = [arb(fmpq(int(i), 1 << pilot['middle_weight_bits']))
               for i in pilot['middle_weight_numerators']]
    if len(weights) != M-N or not all(w > 0 for w in weights):
        raise ValueError('Positive complete archived weights required.')

    # Omitted rows of Z*WZ are PSD. The selected rows alone suffice for
    # the lower bound on K=Z*WZ-F; FLINT verifies their eigenvalue boxes.
    F = arb_mat([[arithmetic.even_entry(i,j,data) for j in selected]
                 for i in selected])
    weighted_F = arb_mat([[weights[i]*F[i,j] for j in range(count)]
                         for i in range(count)])
    partial_K = F.transpose()*weighted_F-F
    eigenvalues = partial_K.eig()
    kappa = 40
    if len(eigenvalues) != count or not all(z.real > kappa for z in eigenvalues):
        raise ArithmeticError('The partial Gram did not certify K>=40I.')
    print('certified K >=', kappa, 'I; seconds', time.monotonic()-started, flush=True)

    # Q(v-Ex) uses the existing complete Fourier assembly, including all
    # pole, Gamma and fourteen prime-power contributions.
    q_u = arb(0)
    last = selected[-1]
    for i in range(last+1):
        q_u += u[i]**2*arithmetic.even_entry(i,i,data)
        q_u += 2*u[i]*sum((u[j]*arithmetic.even_entry(i,j,data)
                          for j in range(i)), arb(0))
    print('finite Q(v-Ex)', q_u.str(24), 'seconds', time.monotonic()-started, flush=True)

    middle_debit = arb(0)
    gradient = arb_mat(count, 1)
    selected_residual = []
    for first in range(N+1, M+1, 128):
        modes = range(first, min(first+128, M+1))
        block = arb_mat([[arithmetic.even_entry(m,n,data)
                          for n in range(last+1)] for m in modes])
        residual = block*ucolumn
        weighted_residual = arb_mat([[weights[m-N-1]*residual[i,0]]
                                    for i,m in enumerate(modes)])
        middle_debit += (residual.transpose()*weighted_residual)[0,0]
        Z = arb_mat([[block[i,N+1+j] for j in range(count)]
                     for i in range(len(modes))])
        gradient += Z.transpose()*weighted_residual
        selected_residual.extend(residual[i,0] for i,m in enumerate(modes) if m <= last)
        if (first-N-1) % 4096 == 0:
            print('middle through', modes.stop-1, 'seconds', time.monotonic()-started, flush=True)
    if len(selected_residual) != count:
        raise ArithmeticError('Incomplete selected residual.')

    # Apply the archived moment core to this single combined vector and
    # its selected gradient columns. No moments or old Gram are regenerated.
    alpha, D, P, H, _ = reconstruct_factors(
        dict(far, N=last), dict(scalar, rows=[rows[n] for n in range(last+1)]))
    from flint import fmpq_mat
    feature = P*fmpq_mat([[x] for x in uq])
    Hfeature = H*feature
    far_center = alpha*D*(feature.transpose()*Hfeature)[0,0]
    bounds_input = SimpleNamespace(
        M=M, L=arithmetic.d, pi=arithmetic.pi,
        budget=2*arithmetic.H+sum((w for w,_ in arithmetic.primes),arb(0)))
    scalar_radius = max(arb(int(rows[n]['s'][1])-int(rows[n]['s'][0])) /
                        (1 << (scalar['grid_bits']+1)) for n in range(last+1))
    bounds = far_bounds(bounds_input, last, far['jet_order'], scalar_radius,
                        far['grid_bits'], far['moment_diagonals'])
    epsilon = arb(D)*bounds['unweighted_matrix_error']
    far_debit = arb(far_center)+arb(0, (epsilon*arb(norm_u_sq)).upper())
    gradient_error = epsilon*arb(norm_u_sq).sqrt()
    for j in range(count):
        center = alpha*D*sum((P[i,N+1+j]*Hfeature[i,0]
                              for i in range(P.nrows())), fmpq(0))
        gradient[j,0] += arb(center)+arb(0, gradient_error.upper())-selected_residual[j]

    fixed_margin = q_u-middle_debit-far_debit
    gradient_sq = sum((gradient[j,0]**2 for j in range(count)), arb(0))
    correction = gradient_sq/kappa
    optimal_upper = fixed_margin+correction
    if not optimal_upper < 0:
        raise ArithmeticError('No strictly negative optimized-family upper bound.')
    result = {
        'schema':f'ACTUAL_C27_K{count}_OBSTRUCTION_V1', 'N':N, 'M':M,
        'selected_exterior_columns':list(selected), 'basis':scalar['basis'],
        'working_bits':args.precision, 'grid_bits':args.output_bits,
        'candidate_sha256':hashlib.sha256(args.candidate.read_bytes()).hexdigest(),
        'candidate_norm_sq':str(norm_v_sq), 'combined_norm_sq':str(norm_u_sq),
        'K_lower_bound':kappa, 'partial_K_rows':list(selected),
        'partial_K_eigenvalue_real_intervals':[grid_interval(z.real,args.output_bits) for z in eigenvalues],
        'finite_combined_Q':grid_interval(q_u,args.output_bits),
        'middle_residual_debit':grid_interval(middle_debit,args.output_bits),
        'far_residual_debit':grid_interval(far_debit,args.output_bits),
        'extended_far_error_grid_bits':640,
        'extended_far_matrix_error':grid_interval(epsilon,640),
        'gradient_intervals':[grid_interval(gradient[j,0],args.output_bits) for j in range(count)],
        'fixed_trial_margin':grid_interval(fixed_margin,args.output_bits),
        'optimizer_correction_upper':grid_interval(correction,args.output_bits),
        'optimized_family_quadratic_upper':grid_interval(optimal_upper,args.output_bits),
        'strictly_negative_upper':True, 'elapsed_seconds':time.monotonic()-started,
        'python_flint':flint.__version__, 'FLINT':flint.__FLINT_VERSION__,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'input_sha256':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['arithmetic-scalars.json.zlib','diagonal-floor-witness.json','middle-scalars.json.zlib',
             'far-coupling.json.zlib','inverse-residual-pilot.json.zlib']},
        'conditions':['Existing complete Fourier/form identification and actual exterior domain.',
                      'The archived positive diagonal W majorizes the same actual exterior inverse.'],
        'scope':f'Directed conditional exclusion of every X for the fixed{count} exterior columns and '
                'archived W: K>=40I and completion of the square give '
                'v*(A-G+T*K^-1T)v<=fixed_trial_margin+||Kx-Tv||²/40<0. '
                'This excludes this sufficient lower-certificate family, not actual Weil positivity, '
                'other exterior spaces/majorants, cofinal positivity or RH. No Lean verification.'}
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print('fixed margin',fixed_margin.str(30),'correction',correction.str(16),
          'optimized upper',optimal_upper.str(30),flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--propose', action='store_true')
    parser.add_argument('--candidate', type=Path, default=HERE/'inverse-residual-k64-candidate.json')
    parser.add_argument('--precision', type=int, default=768)
    parser.add_argument('--output-bits', type=int, default=512)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('Use a fresh output path.')
    ctx.prec = args.precision
    if args.propose:
        propose(args)
    else:
        certify(args)


if __name__ == '__main__':
    main()
