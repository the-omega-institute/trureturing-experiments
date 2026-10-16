"""Append missing within-bin Gamma rate contrasts to archived parent channels.

Use classical finite orthogonal projection, the existing digamma boundary
columns and finite Cauchy moment products. Parent aggregates, scalar rows,
residual generators and signed far moments are consumed without production.
"""
import argparse
import hashlib
from pathlib import Path
import time

from flint import arb, arb_mat, arb_poly, ctx, fmpq

from arb_arithmetic import Arithmetic
from arb_inverse_residual_obstruction import HERE, M, N, inputs
from diagonal_floor_witness import interval_ball
from gamma_boundary_update import boundary_columns, read_json
from gamma_support_optimization import setup
from joint_gamma_boundary import matrix_in, matrix_out, write_packet
from structured_inverse_residual_interval import envelope


def refinement_tree(bins, depth):
    leaves, details = [], []

    def split(root, a, b, remaining, level):
        if b-a == 1 or remaining == 0:
            leaves.append((root, a, b))
            return
        middle = (a+b)//2
        details.append(dict(root=root, left=a, split=middle, right=b, level=level))
        split(root, a, middle, remaining-1, level+1)
        split(root, middle, b, remaining-1, level+1)

    for root, (a, b) in enumerate(bins):
        split(root, a, b, depth, 1)
    return leaves, details


def contrast_columns(arithmetic, bins, depth, debit):
    leaves, details = refinement_tree(bins, depth)
    raw, _ = boundary_columns(arithmetic, [(a, b) for _, a, b in leaves], debit)
    parent_count = len(bins)
    transform = arb_mat(len(leaves), parent_count+len(details))
    c = lambda a: 1-(-arb(fmpq(4*a+1, 2))*arithmetic.d).exp()
    for row, (root, a, b) in enumerate(leaves):
        root_a, root_b = bins[root]
        rescale = (c(root_a)/c(a)).sqrt()
        transform[row, root] = rescale*arb(fmpq(b-a, root_b-root_a)).sqrt()
        for j, detail in enumerate(details):
            if detail['root'] != root:
                continue
            lo, mid, hi = detail['left'], detail['split'], detail['right']
            left, right, total = mid-lo, hi-mid, hi-lo
            if lo <= a and b <= mid:
                coefficient = arb(fmpq((b-a)*right, left*total)).sqrt()
            elif mid <= a and b <= hi:
                coefficient = -arb(fmpq((b-a)*left, right*total)).sqrt()
            else:
                coefficient = arb(0)
            transform[row, parent_count+j] = rescale*coefficient
    combined = raw*transform
    parents = arb_mat([[combined[i, j] for j in range(parent_count)] for i in range(M-N)])
    new = arb_mat([[combined[i, parent_count+j] for j in range(len(details))] for i in range(M-N)])
    # b/(b^2+omega^2) is 1/omega^2-Lipschitz in b. The two
    # arithmetic-progression means differ by the parent interval length.
    constants = []
    for detail in details:
        a, mid, b = detail['left'], detail['split'], detail['right']
        root_a = bins[detail['root']][0]
        constants.append(2*(c(root_a)/arithmetic.d).sqrt()
                         *arb((mid-a)*(b-mid)*(b-a)).sqrt()
                         *(arithmetic.d/(2*arithmetic.pi))**2)
    eta = [debit*x*x/(3*M**3) for x in constants]
    return parents, new, eta, leaves, details


def cauchy_projection(weighted, row_generator, column_generator, modes, columns,
                      kappa, scale, order):
    """Existing finite disjoint-node geometric product, with its absolute tail."""
    count = weighted.ncols()
    moments, absolute = arb_mat(order, 2*count), arb_mat(1, 2*count)
    for first in range(0, len(modes), 256):
        block = modes[first:first+256]
        powers, values = arb_mat(len(block), order), arb_mat(len(block), 2*count)
        for i, m in enumerate(block):
            ratio, power = arb(fmpq(scale*scale, m*m)), arb(1)
            for r in range(order):
                powers[i, r] = power
                power *= ratio
            gm = [row_generator(m, j) for j in range(2)]
            for k in range(count):
                for j in range(2):
                    value = weighted[m-N-1, k]*gm[j]/(m*m)
                    values[i, 2*k+j] = value
                    absolute[0, 2*k+j] += abs(value)
        moments += powers.transpose()*values
    ratio = arb(fmpq(scale, modes[0]))
    remainder = ratio**(2*order)/(1-ratio*ratio)
    result = arb_mat(count, len(columns))
    for k in range(count):
        pa = arb_poly([moments[r, 2*k] for r in range(order)])
        pb = arb_poly([moments[r, 2*k+1] for r in range(order)])
        for j, n in enumerate(columns):
            t = arb(fmpq(n*n, scale*scale))
            error = kappa*remainder*(abs(column_generator[j, 1])*absolute[0, 2*k]
                                     +abs(column_generator[j, 0])*absolute[0, 2*k+1])
            result[k, j] = kappa*(column_generator[j, 1]*pa(t)-column_generator[j, 0]*pb(t))
            result[k, j] += arb(0, error.upper())
    info = dict(order=order, first_mode=modes[0], last_mode=modes[-1],
                ratio=envelope(ratio, 640), remainder=envelope(remainder, 640), grid_bits=640)
    return result, info


def produce(args):
    started = time.monotonic()
    old, _, _, generators, boundary, selected, GE, _, F, _, _, provenance = setup(args)
    parent_action, parent_hash = read_json(args.parent_finite)
    if (parent_action['schema'], parent_action['metadata']) != (
            'ACTUAL_C27_GAMMA_SUPPORT_ACTION_V1', provenance):
        raise ValueError('Require the matching archived parent support action.')
    scalar, rows, pilot, _ = inputs()
    arithmetic = Arithmetic(27, scalar['grid_bits'])
    weights = [arb(fmpq(int(x), 1 << pilot['middle_weight_bits']))
               for x in pilot['middle_weight_numerators']]
    debit = arb(fmpq(int(pilot['far_weight_numerator']), int(pilot['far_weight_denominator'])))
    U, V, eta, leaves, details = contrast_columns(arithmetic, boundary['bins'], args.depth, debit)
    WV = arb_mat([[weights[i]*V[i, k] for k in range(V.ncols())] for i in range(M-N)])
    cross, gram = U.transpose()*WV, V.transpose()*WV
    print('new orthogonal rate-contrast Gram/cross assembled; count', V.ncols(),
          'seconds', time.monotonic()-started, flush=True)
    GO = matrix_in(generators['GO'], generators['grid_bits'])
    GP = arb_mat([[interval_ball(x, old['grid_bits']) for x in row[:2]] for row in old['generators']])
    kappa = interval_ball(old['kappa'], old['grid_bits'])
    last = selected[-1]
    beta, beta_info = cauchy_projection(WV, lambda m, j: GO[m-last-1, j], GP,
        list(range(last+1, M+1)), list(range(N+1)), kappa, N, 512)
    print('new full-coefficient rate-contrast response assembled; seconds',
          time.monotonic()-started, flush=True)
    s = {n: interval_ball(rows[n]['s'], scalar['grid_bits']) for n in range(N+1, M+1)}
    original = lambda m, j: arb((-1)**m)*(m*s[m] if j else 1)
    support, support_info = cauchy_projection(WV, original, GE,
        list(range(2*last+1, M+1)), selected, kappa, last, 384)
    for first in range(N+1, 2*last+1, 128):
        modes = list(range(first, min(first+128, 2*last+1)))
        Z = arb_mat(len(modes), len(selected))
        for i, m in enumerate(modes):
            gm0, gm1 = original(m, 0), original(m, 1)
            for j, n in enumerate(selected):
                Z[i, j] = (F[j, j] if m == n else
                    kappa*(gm0*GE[j, 1]-gm1*GE[j, 0])/(m*m-n*n))
        support += arb_mat([[WV[m-N-1, k] for k in range(V.ncols())] for m in modes]).transpose()*Z
    result = dict(schema='ACTUAL_C27_GAMMA_RATE_CONTRAST_V1', metadata=provenance,
                  parent_support_action_sha256=parent_hash, grid_bits=args.output_bits,
                  working_bits=args.precision, depth=args.depth, leaves=leaves, contrasts=details,
                  parent_channel_count=U.ncols(), new_channel_count=V.ncols(),
                  finite_cross=matrix_out(cross, args.output_bits),
                  finite_gram=matrix_out(gram, args.output_bits),
                  finite_beta=matrix_out(beta, args.output_bits),
                  finite_support=matrix_out(support, args.output_bits),
                  far_gram_trace_bound=envelope(sum(eta, arb(0)), args.output_bits),
                  far_channel_energy_bounds=[envelope(x, args.output_bits) for x in eta],
                  beta_geometry=beta_info, support_geometry=support_info,
                  support_direct_last_mode=2*last,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  boundary_supplier_sha256=hashlib.sha256((HERE/'gamma_boundary_update.py').read_bytes()).hexdigest(),
                  elapsed_seconds=time.monotonic()-started,
                  scope='Only new within-bin orthogonal Gamma rate contrasts and their actual common finite actions. '
                    'Classical Bessel/Haar projection; no new generic theorem or novelty claim. '
                    'Old means, Gram/Beta/support aggregates and scientific scalar/far suppliers are reused. '
                    'No claim of retained positivity, cofinal positivity, Lean verification or RH.')
    write_packet(args.out, result)
    print('new rate-contrast common actions archived; seconds', time.monotonic()-started, flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--target', type=Path, default=HERE/'structured-joint-2048-target.json.zlib')
    p.add_argument('--pivots', type=Path, default=HERE/'structured-joint-2048-pivots.json')
    p.add_argument('--support', type=Path, default=HERE/'support-2048-finite-action.json.zlib')
    p.add_argument('--generators', type=Path, default=HERE/'joint-residual-generators-2048.json.zlib')
    p.add_argument('--boundary', type=Path, default=HERE/'joint-gamma-6765-finite.json.zlib')
    p.add_argument('--parent-finite', type=Path, default=HERE/'gamma-support-6765-finite.json.zlib')
    p.add_argument('--depth', type=int, default=1)
    p.add_argument('--precision', type=int, default=1024)
    p.add_argument('--output-bits', type=int, default=640)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    if args.out.exists() or args.depth < 1 or args.precision < max(640, args.output_bits)+192:
        p.error('Use a fresh output, positive refinement depth and at least192 working guard bits.')
    ctx.prec = args.precision
    produce(args)


if __name__ == '__main__':
    main()
