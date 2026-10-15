"""A certified low-rank input for the actual c27 infinite far coupling.

Uses the existing arithmetic scalars and normalized moment supplier. A finite
geometric expansion is applied to the original divided difference; it keeps
every boundary feature and signed moment. Explicit operator-norm allowances
pay geometric, moment, feature and prefactor errors. Middle modes and the
retained sign are outside this input.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import zlib

from flint import arb, ctx, fmpq, fmpq_mat

from arb_arithmetic import grid_interval
from arb_tail_moments import TailMoments


HERE = Path(__file__).resolve().parent
SCALAR_SHA256 = 'f993be1aa48824cbf22276ac2d80e834ca9e54e4519e2040f437ebb4f724990f'
SEED_SHA256 = 'fc91b487a48f82cb276c0823f28c4a38432f12f822d6e0ed2644a67da7c22ad1'


def pinned_packet(path, expected):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError('Input packet does not match its published hash.')
    return json.loads(zlib.decompress(raw))


def exact_midpoint(pair, bits):
    return fmpq(int(pair[0])+int(pair[1]), 1 << (bits+1))


def upper_dyadic(value, bits):
    if not value.is_finite() or not value >= 0:
        raise ArithmeticError('Nonfinite or negative error bound.')
    return str(int((value.upper()*(1 << bits)).ceil().unique_fmpz()))


def far_bounds(arithmetic, N, R, scalar_radius, bits, moment_diagonals=None):
    """Bounds for the actual infinite Gram, independent of any retained sign."""
    M, B, pi = arithmetic.M, arithmetic.budget, arithmetic.pi
    if not 0 <= N < M or R < 1:
        raise ValueError('Use 0<=N<M and jet order>=1.')
    d = R if moment_diagonals is None else moment_diagonals
    if type(d) is not int or not 1 <= d <= 2*R-1:
        raise ValueError('Use 1<=moment-diagonals<=2*jet-order-1.')
    q, kappa = arb(N)/M, 1-arb(N)/(M+1)
    C2 = B*B/(pi*pi*M)*(4*N/kappa**2+2)
    r2 = 4*B*B*N/(pi*pi*kappa**2*M*(4*R+1))*q**(4*R)
    geometric = 2*(C2*r2).sqrt()+r2

    # Pay for omitted k+l>=d core entries after exact q^(2k) scaling.
    # The sparse core need not be PSD. The scalar operator allowance
    # suffices; no matrix-sign or PSD-completion algorithm is required.
    scaled_P2 = arb(fmpq(1, 2))+R*N*(1+B*B*q*q)
    core_tail = (arb(0) if d == 2*R-1 else
                 q**(2*d)/(1-q*q)*(B*B/(2*d+1)+B/(2*d+2)))
    core_error = 4/(pi*pi*M)*scaled_P2*core_tail

    # The exact feature P has A_0(0)=1/sqrt(2), A_k(n)=(n/M)^(2k),
    # B_k(n)=s_n*(n/M)^(2k+1), with all other index-zero entries zero.
    P2 = arb(fmpq(1, 2))+N*(1+B*B*q*q)/(1-q**4)
    root2 = grid_interval(1/arb(2).sqrt(), bits)
    root2_radius = arb(int(root2[1])-int(root2[0]))/(1 << (bits+1))
    p_error = (root2_radius**2+scalar_radius**2*N*q*q/(1-q**4)).sqrt()
    # H is the genuine PSD moment Gram. Entry masking decreases its
    # Frobenius norm, so the sparse core also has norm <= trace(H).
    Hnorm = sum((B*B/(4*k+1)+arb(1)/(4*k+3) for k in range(R)), arb(0))
    # Every exported moment has width<=2 grid units; all callers use its
    # one exact midpoint, not separately favorable choices of endpoints.
    Herror = 2*R*arb(2)**(-bits)
    p0_norm = P2.sqrt()+p_error
    feature_moment = Hnorm*(2*P2.sqrt()*p_error+p_error**2)+p0_norm**2*Herror
    alpha = 4/(pi*pi*M)
    alpha_pair = grid_interval(alpha, bits)
    alpha_radius = arb(int(alpha_pair[1])-int(alpha_pair[0]))/(1 << (bits+1))
    rounding = alpha*feature_moment+alpha_radius*p0_norm**2*(Hnorm+Herror)
    # Existing c27 joint prime bound59/10 and Neumann-Gamma comparison.
    # Monotonic h_m supplies the constant sufficient inverse weight.
    hA = (arb(M+1)/arithmetic.L).log()-arithmetic.L/(pi*(M+1))-arb(fmpq(59, 10))
    if not hA > 0:
        raise ValueError('The existing c27 high-complement floor is not positive here.')
    debit_grid = 256
    debit_numerator = upper_dyadic(1/hA, debit_grid)
    debit = arb(fmpq(int(debit_numerator), 1 << debit_grid))
    weighted_error = debit*(geometric+core_error+rounding)
    return {
        'coupling_frobenius_sq': C2, 'jet_remainder_frobenius_sq': r2,
        'geometric_gram_error': geometric, 'feature_frobenius_sq': P2,
        'scaled_feature_frobenius_sq': scaled_P2,
        'scaled_moment_tail_norm': core_tail,
        'moment_truncation_gram_error': core_error,
        'feature_frobenius_error': p_error, 'moment_gram_error': Herror,
        'moment_gram_norm': Hnorm, 'rounding_gram_error': rounding,
        'h_first_far': hA, 'unweighted_matrix_error': geometric+core_error+rounding,
        'weighted_matrix_error': weighted_error,
        'alpha_pair': alpha_pair, 'zero_feature_pair': root2,
        'debit_numerator': debit_numerator, 'debit_grid_bits': debit_grid}


def reconstruct_factors(packet, scalar, centered=True):
    """Exact rational alpha,D,P0,H0,error for the published far majorant.

    The centered unphased convention multiplies retained column n by (-1)^n.
    No dense (N+1)-square matrix or matrix-sign algorithm is constructed.
    """
    N, M, R, bits = packet['N'], packet['M'], packet['jet_order'], packet['grid_bits']
    d = packet['moment_diagonals']
    if packet['schema'] != 'ACTUAL_C27_FAR_COUPLING_GRAM_V1' or scalar['c'] != 27:
        raise ValueError('Use the declared actual c27 inputs.')
    low = {row['n']: exact_midpoint(row['s'], scalar['grid_bits']) for row in scalar['rows']}
    P = [[fmpq(0) for _ in range(N+1)] for _ in range(2*R)]
    P[0][0] = exact_midpoint(packet['zero_feature_pair'], bits)
    for n in range(1, N+1):
        x, phase = fmpq(n, M), (-1 if centered and n % 2 else 1)
        for k in range(R):
            P[k][n] = phase*x**(2*k)
            P[R+k][n] = phase*low[n]*x**(2*k+1)
    U = {row['U_power']: exact_midpoint(row['U'], bits) for row in packet['moment_pairs']}
    T = {row['T_power']: exact_midpoint(row['T'], bits) for row in packet['moment_pairs']}
    Z = {row['p']: exact_midpoint(row['interval'], bits) for row in packet['zeta']}
    H = [[fmpq(0) for _ in range(2*R)] for _ in range(2*R)]
    for k in range(R):
        for l in range(R):
            if k+l >= d:
                continue
            H[k][l] = U[2*k+2*l+2]
            H[k][R+l] = H[R+l][k] = -T[2*k+2*l+3]
            H[R+k][R+l] = Z[2*k+2*l+4]
    return (exact_midpoint(packet['alpha_pair'], bits),
            fmpq(int(packet['D_numerator']), 1 << packet['D_grid_bits']),
            fmpq_mat(P), fmpq_mat(H),
            fmpq(int(packet['upper_bounds']['weighted_matrix_error']),
                 1 << packet['bound_grid_bits']))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scalar-input', type=Path, default=HERE/'arithmetic-scalars.json.zlib')
    parser.add_argument('--seed-input', type=Path, default=HERE/'tail-moments.json.zlib')
    parser.add_argument('--last-low', type=int, default=1280)
    parser.add_argument('--cutoff', type=int, default=32768)
    parser.add_argument('--jet-order', type=int, default=56)
    parser.add_argument('--stirling-order', type=int, default=32)
    parser.add_argument('--precision', type=int, default=2304)
    parser.add_argument('--full-moments', action='store_true',
                        help='Supply all core diagonals; highest powers need more precision.')
    parser.add_argument('--output-bits', type=int, default=640)
    parser.add_argument('--matrix-bits', type=int, default=500)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--moment-cache', type=Path,
                        help='Save completed moment pairs atomically for reuse.')
    parser.add_argument('--resume-sha256',
                        help='Required SHA-256 when reusing an existing moment cache.')
    args = parser.parse_args()
    if args.out.exists():
        parser.error('Use a fresh output path.')
    if args.moment_cache and args.moment_cache.resolve() == args.out.resolve():
        parser.error('Moment cache and final output must be different paths.')
    if args.output_bits < 64 or args.precision < args.output_bits+64 or args.matrix_bits < 1:
        parser.error('Use output-bits>=64, guard bits>=64 and matrix-bits>=1.')
    ctx.prec = args.precision
    scalar = pinned_packet(args.scalar_input, SCALAR_SHA256)
    seed = pinned_packet(args.seed_input, SEED_SHA256)
    N, M, R = args.last_low, args.cutoff, args.jet_order
    d = 2*R-1 if args.full_moments else R
    if scalar['c'] != 27 or not 0 <= N <= scalar['last_low'] or not N < M or R < 1:
        parser.error('Require c27 scalar coverage, 0<=N<M and jet-order>=1.')
    low = {row['n']: row for row in scalar['rows'] if row['n'] <= N}
    if set(low) != set(range(N+1)):
        raise ValueError('Missing actual low symbol input.')
    scalar_radius = arb(max(int(row['s'][1])-int(row['s'][0]) for row in low.values()))
    scalar_radius /= 1 << (scalar['grid_bits']+1)
    arithmetic = TailMoments(27, M, args.stirling_order)
    bounds = far_bounds(arithmetic, N, R, scalar_radius, args.output_bits, d)
    if not bounds['weighted_matrix_error'] < arb(2)**(-args.matrix_bits):
        raise ArithmeticError('Jet/rounding allowance exceeds requested matrix accuracy.')
    print('matrix-error log2 upper',
          (bounds['weighted_matrix_error'].log()/arb(2).log()).upper().str(16), flush=True)

    seed_rows = {row['p']: row for row in seed['rows']}
    use_seed = (seed['c'] == 27 and seed['M'] == M
                and seed['grid_bits'] == args.output_bits)
    started, rows = time.monotonic(), []
    source_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    supplier_hash = hashlib.sha256((HERE/'arb_tail_moments.py').read_bytes()).hexdigest()
    cache_header = {
        'schema': 'ACTUAL_C27_FAR_MOMENT_CACHE_V1', 'c': 27, 'M': M,
        'grid_bits': args.output_bits, 'working_bits': args.precision,
        'stirling_order': args.stirling_order,
        'source_sha256': source_hash, 'moment_supplier_sha256': supplier_hash,
        'seed_packet_sha256': SEED_SHA256}
    cached = {}
    if args.moment_cache and args.moment_cache.exists():
        if not args.resume_sha256:
            parser.error('Pin the existing moment cache with --resume-sha256.')
        saved = pinned_packet(args.moment_cache, args.resume_sha256)
        if saved['header'] != cache_header:
            raise ValueError('Moment cache belongs to different sources or parameters.')
        for row in saved['moment_pairs']:
            p = row['U_power']
            if (set(row) != {'U_power', 'T_power', 'U', 'T'} or
                    type(p) is not int or p < 2 or p % 2 or row['T_power'] != p+1
                    or p in cached):
                raise ValueError('Invalid or duplicate cached moment powers.')
            for family in ('U', 'T'):
                pair = row[family]
                if (not isinstance(pair, list) or len(pair) != 2 or
                        any(type(v) is not str or str(int(v)) != v for v in pair) or
                        not 0 <= int(pair[1])-int(pair[0]) <= 2):
                    raise ValueError('Invalid cached dyadic moment interval.')
            cached[p] = row
    elif args.resume_sha256:
        parser.error('--resume-sha256 requires an existing --moment-cache.')
    reused = 0
    # Check the highest powers first, before paying for all lower moments.
    order = [d-1]+list(range(d-1))
    for i in order:
        up, tp = 2*i+2, 2*i+3
        if up in cached:
            rows.append(cached[up])
            reused += int(use_seed and up in seed_rows)+int(use_seed and tp in seed_rows)
            print('reused cached powers', up, tp, flush=True)
            continue
        row = {'U_power': up, 'T_power': tp}
        for family, p in [('U', up), ('T', tp)]:
            if use_seed and p in seed_rows:
                row[family] = seed_rows[p][family]
                reused += 1
            else:
                value, _ = (arithmetic.square_moment(p) if family == 'U'
                            else arithmetic.linear_moment(p))
                row[family] = grid_interval(value, args.output_bits)
        rows.append(row)
        cached[up] = row
        if args.moment_cache:
            cache_raw = (json.dumps({'header': cache_header,
                                    'moment_pairs': sorted(cached.values(),
                                                           key=lambda r: r['U_power'])},
                                   sort_keys=True, separators=(',', ':'))+'\n').encode()
            temporary = args.moment_cache.with_name(args.moment_cache.name+'.tmp')
            temporary.write_bytes(zlib.compress(cache_raw, level=9))
            temporary.replace(args.moment_cache)
        print('completed powers', up, tp, 'kernels', len(arithmetic._kernels), flush=True)
    rows.sort(key=lambda r: r['U_power'])
    # No new zeta algorithm: the same existing Hurwitz supplier is used.
    zeta = [{'p': 2*i+4, 'interval': grid_interval(arithmetic.kernel(2*i+4).real,
                                                args.output_bits)}
            for i in range(d)]
    bound_bits = args.matrix_bits+128
    result = {
        'schema': 'ACTUAL_C27_FAR_COUPLING_GRAM_V1',
        'c': 27, 'N': N, 'M': M, 'first_far': M+1, 'jet_order': R,
        'moment_diagonals': d, 'stirling_order': args.stirling_order,
        'feature_dimension': 2*R, 'working_bits': args.precision,
        'grid_bits': args.output_bits, 'bound_grid_bits': bound_bits,
        'matrix_error_bits': args.matrix_bits,
        'scalar_packet_sha256': SCALAR_SHA256, 'seed_packet_sha256': SEED_SHA256,
        'reused_scalar_moments': reused, 'moment_pairs': rows, 'zeta': zeta,
        'feature_definition': {
            'A_k(0)': 'midpoint(zero_feature_pair) if k=0; otherwise0',
            'A_k(n)': '(n/M)^(2*k), n>=1',
            'B_k(n)': 'midpoint(s_n)*(n/M)^(2*k+1), n>=1; zero at n=0',
            'scalar_midpoint_denominator_bits': scalar['grid_bits']+1},
        'moment_gram_definition': {
            'AA_kl': 'midpoint(U_(2*k+2*l+2))',
            'AB_kl': '-midpoint(T_(2*k+2*l+3))',
            'BB_kl': 'midpoint(Z_(2*k+2*l+4))',
            'mask': 'All blocks are zero where k+l>=moment_diagonals.'},
        'unweighted_majorant': 'midpoint(alpha_pair)*P0^T*H0*P0 + error*I',
        'weighted_majorant': 'D*(midpoint(alpha_pair)*P0^T*H0*P0) + weighted_error*I',
        'unphased_basis_transport': 'Multiply column n of P0 by (-1)^n.',
        'alpha_pair': bounds.pop('alpha_pair'), 'zero_feature_pair': bounds.pop('zero_feature_pair'),
        'h_first_far_interval': grid_interval(bounds.pop('h_first_far'), bound_bits),
        'D_numerator': bounds.pop('debit_numerator'), 'D_grid_bits': bounds.pop('debit_grid_bits'),
        'upper_bounds': {k: upper_dyadic(v, bound_bits) for k, v in bounds.items()},
        'kernel_count': len(arithmetic._kernels),
        'source_sha256': source_hash,
        'moment_supplier_sha256': supplier_hash,
        'elapsed_seconds': time.monotonic()-started,
        'scope': 'Actual c27 infinite far coupling only, conditional on existing '
                 'Fourier/form and Neumann-Gamma identifications and prime bound59/10. '
                 'No middle modes, retained sign, complete coupling, Lean or RH certificate.'}
    if int(result['upper_bounds']['weighted_matrix_error']) >= 1 << (bound_bits-args.matrix_bits):
        raise ArithmeticError('Outward exported allowance exceeds matrix target.')
    raw = (json.dumps(result, sort_keys=True, separators=(',', ':'))+'\n').encode()
    args.out.write_bytes(zlib.compress(raw, level=9))
    print('written', str(args.out), 'packet_sha256', hashlib.sha256(args.out.read_bytes()).hexdigest())


if __name__ == '__main__':
    main()
