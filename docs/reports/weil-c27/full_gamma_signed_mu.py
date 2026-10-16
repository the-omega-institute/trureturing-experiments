"""Enclose the actual same-source coherent signed proxy functional.

Reuse the saved full56 short jet and all original source certificates.
FLINT supplies reciprocal-sum convolution, polynomial evaluation and
real Hurwitz kernels. Abel summation bounds the oscillatory infinite tail.
Only the source-specific assembly and errors are new.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time
import zlib

from flint import arb, arb_poly, ctx, fmpq

from arb_tail_moments import TailMoments, phase_representative
from diagonal_floor_witness import interval_ball
from full_gamma_coherent_split_capacity import freeze_upper
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS, M, FIRST, SPLIT = 640, 32768, 1281, 2048
MAX_ORDER = 256
GEOMETRIC_BUDGET = fmpq(1, 64)
PINS = {
    'full-gamma-signed-jet-interface.json.zlib':
        '154949d8ab3e1f2875de531aa63a9c5a3c3d7d133640a5ce8b1d23c98d67565a',
    'tail-moments.json.zlib':
        'fc91b487a48f82cb276c0823f28c4a38432f12f822d6e0ed2644a67da7c22ad1',
    'full-gamma-dual-32768-finite.json.zlib':
        '4085fa7048bafc987c75b2340b56d98d556534ae74d304e61469d4e6e52af1b5',
    'full-gamma-dual-32768-defect.json.zlib':
        '6f320576338562a4e7338906ca66b25ab49e4a1d4460e87e26165cd9579a1154',
    'full-gamma-primal-tail.json':
        'fad634fb25d8d20e83a15d8277569978646246c2a560d0dc0f7a854b5087455a',
    'full-gamma-hankel-tail-energy.json.zlib':
        'd09a5ab254d164cfd40982b31e1a155e31d6288b705d916bf248236a0d527606',
}


def require_finite(values):
    if any(not value.is_finite() for value in values):
        raise ArithmeticError('Nonfinite signed-functional value; stop without retry.')


def cached_symbol(seed, scale, first_tail):
    """Load published coefficients; initialize only needed phases/kernels."""
    supplier = object.__new__(TailMoments)
    supplier.c, supplier.M, supplier.A = seed['c'], scale, first_tail
    supplier.pi, supplier.L = arb.pi(), arb(seed['c']).log()
    supplier._kernels, supplier._phases = {}, {}
    supplier.primes = {}
    for r, p in seed['prime_powers']:
        u, sign = phase_representative(Fraction(r), supplier.c)
        supplier.primes[u] = supplier.primes.get(u, arb(0))+sign*arb(p).log()/arb(r).sqrt()
    supplier.primes = {u: w for u, w in supplier.primes.items() if not (w == 0)}
    if len(supplier.primes) != seed['folded_sine_phases']:
        raise ValueError('Different folded actual-symbol phases.')
    coefficients = [interval_ball(pair, seed['grid_bits'])
                    for pair in seed['normalized_polynomial']]
    sigma = interval_ball(seed['uniform_symbol_error'], seed['analytic_grid_bits']).upper()
    require_finite(coefficients+[sigma]+list(supplier.primes.values()))
    if sigma < 0 or len(coefficients) != 64 or not coefficients[0].overlaps(-arb.pi()/4):
        raise ValueError('Different certified smooth polynomial or symbol error.')
    return supplier, arb_poly(coefficients), sigma


def reciprocal_sum_block(z, first_source, first_output, count):
    """One existing FLINT product for disjoint reciprocal-sum samples."""
    if not z or first_source < 1 or first_output <= first_source+len(z)-1 or count < 1:
        raise ValueError('Require a positive disjoint source/output block.')
    require_finite(z)
    size = len(z)
    kernel = arb_poly([arb(fmpq(1, first_output+first_source+a))
                       for a in range(size+count-1)])
    product = arb_poly(list(reversed(z)))*kernel
    result = [product[size-1+i] for i in range(count)]
    require_finite(result)
    return result


def source_polynomial(z, first_source, scale, order):
    """Finite signed power sums for the separated reciprocal remainder."""
    if not z or first_source < 1 or scale < first_source+len(z)-1 or order < 1:
        raise ValueError('Invalid normalized source-power domain.')
    ratios = [arb(fmpq(-(first_source+i), scale)) for i in range(len(z))]
    terms = list(z)
    coefficients = []
    for _ in range(order):
        coefficients.append(sum(terms, arb(0)))
        terms = [v*x for v, x in zip(terms, ratios)]
    require_finite(coefficients)
    return arb_poly(coefficients)


def finite_functional(z, first_source, scale, apoly, bpoly, smooth, supplier, kappa, weight):
    first, count = scale+1, scale
    samples = reciprocal_sum_block(z, first_source, first, count)
    xs = [arb(fmpq(scale, first+i)) for i in range(count)]
    squares = [x*x for x in xs]
    av = apoly.evaluate(squares, algorithm='iter')
    bv = bpoly.evaluate(squares, algorithm='iter')
    symbols = smooth.evaluate(xs, algorithm='iter')
    for u, w in supplier.primes.items():
        phase = supplier.phase(u)
        power = phase**first
        for i in range(count):
            symbols[i] -= w*power.imag
            power *= phase
    require_finite(list(av)+list(bv)+list(symbols))
    brackets = [x*b-s*a for x, a, b, s in zip(xs, av, bv, symbols)]
    total = sum((s*bracket/(first+i) for i, (s, bracket)
                 in enumerate(zip(samples, brackets))), arb(0))*kappa*weight/arb.pi()
    jet = [2*(-1)**(first+i)*bracket/(arb.pi()*(first+i))
           for i, bracket in enumerate(brackets)]
    require_finite([total]+jet)
    return total, samples, jet


def contracted_tail(source, apoly, bpoly, smooth, supplier, kappa, weight):
    # Convert A(x^2),B(x^2) using the existing polynomial representation.
    ac = [arb(0)]*(2*len(apoly)-1)
    bc = [arb(0)]*(2*len(bpoly))
    for k, v in enumerate(apoly.coeffs()):
        ac[2*k] = v
    for k, v in enumerate(bpoly.coeffs()):
        bc[2*k+1] = v
    even_a = arb_poly(ac)
    f0 = source*(arb_poly(bc)-smooth*even_a)
    f1 = source*even_a
    # Arb comparisons are definite relations; only discard provably exact zeros.
    total = sum((v*supplier.kernel(p+2).real
                 for p, v in enumerate(f0.coeffs()) if not (v == 0)), arb(0))
    total *= kappa*weight/(arb.pi()*supplier.M)
    # For the exact unit phase xi != 1, Abel summation gives
    # |sum_{m>=T} xi^m/m^q| <= 2/(|1-xi|*T^q), q>=2.
    # Bound the SAME approximate-source F1; the old telescope is unchanged.
    ratio = arb(fmpq(supplier.M, supplier.A))
    coefficient_cap = sum((abs(v).upper()*ratio**p
                           for p, v in enumerate(f1.coeffs())), arb(0)).upper()
    phase_gaps, phase_cap = {}, arb(0)
    for u, w in supplier.primes.items():
        gap = abs(1-supplier.phase(u)).lower()
        if not gap.is_finite() or not gap > 0:
            raise ArithmeticError('Uncertified nontrivial unit-phase gap; stop without retry.')
        phase_gaps[str(u)] = gap
        phase_cap += abs(w).upper()/gap
    oscillatory_error = freeze_upper(2*kappa*weight*coefficient_cap*phase_cap
                                     /(arb.pi()*supplier.A**2))
    require_finite([total, coefficient_cap, phase_cap, arb(oscillatory_error)])
    return total, oscillatory_error, coefficient_cap, phase_cap, phase_gaps


def verify_guards(guards):
    for path, digest in guards.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise RuntimeError('Read dependency changed: '+str(path))


def write_packet(out, packet):
    raw = zlib.compress(json.dumps(packet, separators=(',', ':'), allow_nan=False).encode(), 9)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('xb') as stream:
        stream.write(raw)
    return hashlib.sha256(raw).hexdigest()


def error_contract(z_upper, jcap, principal, h_a, sigma, kappa, weight, scale):
    values = [arb(v) for v in (z_upper, jcap, principal, sigma, kappa, weight)]
    require_finite(values+list(h_a))
    if (not z_upper >= 0 or not jcap >= 0 or not principal >= 0
            or not sigma >= 0 or not kappa > 0 or not weight > 0 or scale < 2):
        raise ValueError('Invalid signed-functional norm/error domain.')
    constant = freeze_upper(arb(kappa)*arb(jcap)*arb(weight).sqrt()*arb(z_upper)
                            /(2*arb(2*scale).sqrt()*arb(GEOMETRIC_BUDGET)))
    # This exact bit-length choice satisfies 2**order >= constant without a scan.
    order = max(1, int(constant.p).bit_length()-int(constant.q).bit_length()+1)
    if order > MAX_ORDER:
        raise ArithmeticError('Fixed source-power budget exceeded; stop without retry.')
    delta = freeze_upper(arb(weight).sqrt()*arb(z_upper)*arb(2)**(-order)
                         /(2*arb((2*order+1)*2*scale).sqrt()))
    ratio = arb(fmpq(scale, scale+1))
    aabs = sum((abs(v).upper()*ratio**(2*k) for k, v in enumerate(h_a)), arb(0)).upper()
    beta = freeze_upper(2*arb(weight).sqrt()*sigma*aabs/(arb.pi()*arb(scale).sqrt()))
    remainder = freeze_upper(arb(kappa)*(arb(jcap)*arb(delta)
                                         +(arb(principal)+arb(delta))*arb(beta)))
    require_finite([arb(delta), arb(beta), arb(remainder)])
    if not arb(kappa)*arb(jcap)*arb(delta) <= arb(GEOMETRIC_BUDGET):
        raise ArithmeticError('Prescribed geometric error budget not covered.')
    return order, delta, beta, remainder


def consume(directory, out):
    started = time.monotonic()
    source = Path(__file__).resolve()
    invoked = {source}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename and Path(filename).resolve().parent == HERE:
            invoked.add(Path(filename).resolve())
    guards = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in invoked}
    packets = {}
    for name, digest in PINS.items():
        path = directory/name
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted signed-functional packet: '+name)
        guards[path] = digest
        packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    jet, seed, finite, defect, primal, hankel = [packets[name] for name in PINS]
    if [p['schema'] for p in (jet, seed, finite, defect, primal, hankel)] != [
            'ACTUAL_C27_SIGNED_JET_INTERFACE_V1', 'ACTUAL_WEIL_NORMALIZED_TAIL_MOMENTS_V1',
            'ACTUAL_C27_FULL_GAMMA_DUAL_ACTION_V1', 'ACTUAL_C27_FINITE_DUAL_DEFECT_V1',
            'ACTUAL_C27_FULL_GAMMA_PRIMAL_TAIL_V1', 'ACTUAL_C27_SIGNED_HANKEL_TAIL_ENERGY_V1']:
        raise ValueError('Different signed-functional packet roles.')
    if (any(p['grid_bits'] != BITS for p in (jet, seed, finite, defect, primal, hankel))
            or (jet['physical_last_mode'], jet['first_far_mode'], jet['jet_order'])
                != (SPLIT, M+1, 56)
            or jet['feature_order'] != 'A0..A55,B0..B55'
            or (seed['c'], seed['M'], seed['first_mode'], seed['stirling_order'],
                seed['analytic_grid_bits']) != (27, M, M+1, 32, 832)
            or finite['metadata'] != defect['metadata'] or finite['metadata'] != primal['metadata']
            or finite['metadata']['input_sha256'] != jet['input_sha256']
            or (finite['metadata']['N'], finite['metadata']['M'],
                finite['metadata']['last_mode'], finite['dual_first'], finite['dual_last'],
                finite['dual_bits']) != (FIRST-1, M, SPLIT, SPLIT+1, M, 384)
            or finite['metadata']['basis'] !=
                '1/sqrt(d),sqrt(2/d)*cos(2*pi*n*x/d); no (-1)^n phase'
            or (jet['source_bits'], jet['source_row']) != (384, 0)
            or any(jet[k] != hankel[k] for k in
                   ('source_bits', 'source_row', 'source_integer_sha256'))
            or (hankel['first_source_mode'], hankel['split_source_mode'],
                hankel['last_source_mode'], hankel['first_far_mode']) != (FIRST, SPLIT, M, M+1)
            or primal['tail_seed_sha256'] != PINS['tail-moments.json.zlib']
            or defect['finite_action_sha256'] != PINS['full-gamma-dual-32768-finite.json.zlib']
            or any(primal[k] != PINS[n] for k, n in [
                ('finite_action_sha256', 'full-gamma-dual-32768-finite.json.zlib'),
                ('finite_defect_sha256', 'full-gamma-dual-32768-defect.json.zlib')])):
        raise ValueError('Different full support, source, symbol, normalization or grids.')
    inherited = dict(jet['inherited_dependency_sha256'])
    for mapping in (jet['invocation_source_sha256'], jet['packet_sha256']):
        for name, digest in mapping.items():
            if name in inherited and inherited[name] != digest:
                raise ValueError('Conflicting common-jet dependency: '+name)
            inherited[name] = digest
    for name, digest in inherited.items():
        path = directory/name
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError('Different inherited producer/helper/input: '+name)
        if path in guards and guards[path] != digest:
            raise ValueError('Different invoked dependency: '+name)
        guards[path] = digest
    # The current TailMoments methods are guarded by the accepted far packet.
    # The seed's historical source SHA is provenance, not a current-source pin.
    qrows, integers = defect['D_inverse_difference'], finite['dual_integers']
    if (len(qrows) != SPLIT-FIRST+1 or any(len(row) != 1 for row in qrows)
            or len(integers) != M-SPLIT
            or len(jet['short_jet_coefficients']) != 112
            or any(len(row) != 1 for row in jet['short_jet_coefficients'])):
        raise ValueError('Incomplete original source or complete short jet.')
    identifier = hashlib.sha256(('\n'.join(str(int(x)) for x in integers)+'\n').encode()).hexdigest()
    if identifier != jet['source_integer_sha256']:
        raise ValueError('Different exact384-bit source integers.')
    exact = {k: fmpq(*map(int, jet['exact_rationals'][k])) for k in
             ('kappa', 'frozen_linear_envelope_upper', 'selected_direction_norm_upper',
              'jet_forward_error_upper', 'signed_proxy_functional_error_upper',
              'strict_mu_upper_threshold')}
    kappa, clin = exact['kappa'], exact['frozen_linear_envelope_upper']
    theta_error = exact['signed_proxy_functional_error_upper']
    weight = fmpq(*map(int, hankel['exact_rationals']['original_far_weight']))
    if (kappa != fmpq(33, 32) or kappa != fmpq(*map(int, hankel['exact_rationals']['kappa']))
            or exact['strict_mu_upper_threshold'] != -clin/2-theta_error
            or not weight > 0 or not theta_error >= 0):
        raise ValueError('Different original weight or exact signed criterion.')
    principal = freeze_upper(interval_ball(
        hankel['scalar_bounds']['principal_weighted_norm_upper'], BITS).upper())
    jcap = exact['selected_direction_norm_upper']+exact['jet_forward_error_upper']
    # gp_*_hankel_term = ||source||_1/(2sqrt(M)), unweighted, in the old producer.
    source_sum = sum((interval_ball(primal['bounds'][key], BITS).upper()
                      for key in ('gp_w_hankel_term', 'gp_q_hankel_term')), arb(0))
    z_upper = freeze_upper(2*arb(M).sqrt()*source_sum)
    h = [interval_ball(row[0], BITS) for row in jet['short_jet_coefficients']]
    z = [-(-1)**(FIRST+i)*interval_ball(row[0], BITS) for i, row in enumerate(qrows)]
    z += [(-1)**(SPLIT+1+i)*arb(fmpq(int(x), 1 << 384)) for i, x in enumerate(integers)]
    require_finite(z+h)
    supplier, smooth, sigma = cached_symbol(seed, M, 2*M+1)
    order, delta, beta, remainder = error_contract(
        z_upper, jcap, principal, h[:56], sigma, kappa, weight, M)
    apoly, bpoly = arb_poly(h[:56]), arb_poly(h[56:])
    print('fixed source-power order', order, 'real kernel limit', order+174,
          'oscillatory Lerch calls', 0, flush=True)
    finite_mu, samples, finite_jet = finite_functional(
        z, FIRST, M, apoly, bpoly, smooth, supplier, kappa, weight)
    print('finite signed functional enclosed; seconds', time.monotonic()-started, flush=True)
    print('finite_approximate_functional', finite_mu.str(24), flush=True)
    cauchy_export = [envelope(v, BITS) for v in samples]
    jet_export = [envelope(v, BITS) for v in finite_jet]
    if any(pair is None or int(pair[0]) > int(pair[1]) for pair in cauchy_export+jet_export):
        raise ArithmeticError('Invalid finite-functional export; stop without retry.')
    checkpoint_path = out.with_suffix('.finite.zlib')
    checkpoint = dict(schema='ACTUAL_C27_SIGNED_MU_FINITE_BLOCK_V1', packet_sha256=PINS,
        source_sha256=guards[source], source_integer_sha256=identifier,
        source_bits=384, source_row=0, grid_bits=BITS, working_bits=1024,
        first_source_mode=FIRST, last_source_mode=M, physical_residual_last_mode=SPLIT,
        finite_first_mode=M+1, finite_last_mode=2*M, normalization_scale=M,
        jet_order=56, feature_order=jet['feature_order'],
        original_far_weight=[str(weight.p), str(weight.q)], kappa=[str(kappa.p), str(kappa.q)],
        finite_approximate_functional=envelope(finite_mu, BITS),
        finite_cauchy_source_outputs=cauchy_export, finite_approximate_jet_outputs=jet_export,
        finite_output_basis='Cauchy=sum z_n/(m+n); approximate jet retains (-1)^m and exact2/pi',
        symbol_is_approximate=True, errors_paid_in_final_packet=True,
        elapsed_seconds=time.monotonic()-started)
    verify_guards(guards)
    checkpoint_sha = write_packet(checkpoint_path, checkpoint)
    print('finite checkpoint saved', checkpoint_path.name, checkpoint_sha, flush=True)
    source_poly = source_polynomial(z, FIRST, M, order)
    tail_mu, osc_error, coefficient_cap, phase_cap, phase_gaps = contracted_tail(
        source_poly, apoly, bpoly, smooth, supplier, kappa, weight)
    alg = finite_mu+tail_mu
    total_error = freeze_upper(arb(remainder)+arb(osc_error))
    mu = alg+arb(0, arb(total_error))
    theta = mu+arb(0, arb(theta_error))
    center = arb(clin)+2*theta
    scalars = dict(finite_approximate_functional=finite_mu, tail_approximate_functional=tail_mu,
        algorithmic_functional=alg, actual_proxy_functional=mu, actual_signed_theta=theta,
        center_upper_expression_before_good_metric=center, cached_source_l1_upper=arb(z_upper),
        source_tail_norm_error_upper=arb(delta), symbol_jet_norm_error_upper=arb(beta),
        functional_analytic_error_upper=arb(remainder), inherited_theta_error_upper=arb(theta_error),
        oscillatory_tail_error_upper=arb(osc_error), complete_functional_error_upper=arb(total_error),
        oscillatory_coefficient_cap=coefficient_cap, inverse_phase_gap_weight_cap=phase_cap,
        unchanged_strict_mu_threshold=arb(exact['strict_mu_upper_threshold']))
    require_finite(list(scalars.values()))
    exports = {}
    for key, value in scalars.items():
        print(key, value.str(24), flush=True)
        exports[key] = envelope(value, BITS)
    if any(pair is None or int(pair[0]) > int(pair[1])
           for pair in list(exports.values())+cauchy_export+jet_export):
        raise ArithmeticError('Invalid actual-functional export; stop without retry.')
    upper_negative = bool(center.upper() < 0)
    verdict = ('negative-upper-before-good-block-metric' if upper_negative else
               'fixed-signed-upper-expression-positive' if center.lower() > 0 else
               'inconclusive-signed-center-upper-expression')
    rationals = dict(original_far_weight=weight, kappa=kappa,
        source_l1_upper=z_upper, jet_norm_upper=jcap, principal_norm_upper=principal,
        geometric_scalar_budget=GEOMETRIC_BUDGET, source_tail_norm_error_upper=delta,
        symbol_jet_norm_error_upper=beta, functional_analytic_error_upper=remainder,
        oscillatory_tail_error_upper=osc_error, complete_functional_error_upper=total_error,
        inherited_theta_error_upper=theta_error, frozen_linear_envelope_upper=clin,
        strict_mu_upper_threshold=exact['strict_mu_upper_threshold'])
    result = dict(schema='ACTUAL_C27_SAME_SOURCE_SIGNED_MU_ABEL_V1', packet_sha256=PINS,
        source_sha256=guards[source], invocation_source_sha256={
            p.name: sha for p, sha in guards.items() if p in invoked},
        inherited_dependency_sha256={p.name: sha for p, sha in guards.items() if p not in invoked},
        source_integer_sha256=identifier, source_bits=384, source_row=0,
        grid_bits=BITS, working_bits=1024, symbol_metadata_grid_bits=seed['analytic_grid_bits'],
        first_source_mode=FIRST, last_source_mode=M, physical_residual_last_mode=SPLIT,
        finite_first_mode=M+1, finite_last_mode=2*M, analytic_first_mode=2*M+1,
        normalization_scale=M, jet_order=56, source_power_order=order,
        feature_order=jet['feature_order'], scalar_bounds=exports,
        exact_rationals={k: [str(v.p), str(v.q)] for k, v in rationals.items()},
        finite_checkpoint=dict(path=checkpoint_path.name, sha256=checkpoint_sha),
        real_hurwitz_kernels={str(p): envelope(v.real, BITS)
                              for (p, u), v in supplier._kernels.items() if u == 1},
        certified_phase_gap_lower={u: envelope(gap, BITS) for u, gap in phase_gaps.items()},
        new_reciprocal_sum_convolutions=1, new_source_power_polynomials=1,
        new_scalar_tail_kernels=len(supplier._kernels), new_phase_initializations=len(supplier._phases),
        new_oscillatory_lerch_kernels=0,
        new_cached_symbol_polynomial_productions=0, new_source_selections=0,
        new_reference_gram_trace_inverse_solve_calls=0, old_producers_replayed=0,
        parameter_source_precision_cutoff_scan=False, actual_same_source_mu_supplied=True,
        actual_signed_theta_interval_supplied=True, negative_upper_before_good_metric=upper_negative,
        good_block_metric_paid=False, actual_full_border_sign_claim=False,
        verdict=verdict, elapsed_seconds=time.monotonic()-started,
        conditions=jet['conditions']+[
            'Finite Cauchy samples use ONE FLINT reciprocal-sum convolution, all q/w/h intervals and the SAME original source. Parities cancel in mu and the scalar prefactor is kappa*D/pi.',
            'The source1-norm cap is reused from the saved UNWEIGHTED gp_w/q_hankel_term=||w/q||1/(2sqrtM), not a weighted norm or energy lower bound.',
            'Only the tail m>=2M+1 uses the finite source expansion. Its norm error includes the principal factor1/2; actual jet norm <=Ncap+eps_b.',
            'Cached smooth polynomial and832-grid uniform symbol error are reused. Only14 elementary prime-power weights,13 folded phases and real Hurwitz endpoint kernels are initialized; no old polynomial or scalar supplier is replayed.',
            'All real kernels have q>=2 at endpoint65537 with normalizationM32768. No complex infinite phase kernel is evaluated. Every exact nontrivial unit phase has a certified positive gap.',
            'For approximate-source F1, Abel summation gives Rosc=up640[2*kappa*D/(pi*T^2)*sum_p|F1_p|(M/T)^p*sum_u|a_u|/|1-xi_u|]. This bounds the omitted oscillatory tail of the SAME approximate functional.',
            'R_mu=kappa*(Jcap*delta_r+(P+delta_r)*beta_s) still pays source-tail truncation and the global symbol error. I_alg +/-(R_mu+Rosc) encloses actual mu; inheritedE_mu additionally encloses actual theta.',
            'Saved finite jet outputs use the approximate symbol, not the actual far operator. Their complete symbol/forward errors remain paid only by the stated norm/error contract.',
            'Clin+2theta is the existing sufficient upper expression before the additional PSD good-block metric. A positive expression does not establish a positive actual center; a negative upper alone does not supply a full-border sign.'
        ], scope='Source-specific actual signed functional using existing FLINT algorithms and inherited certificates. No new generic algorithm or Lean theorem. Good-block metric, actual full-border sign, retained/cofinal positivity, Robin tail and RH remain unresolved.')
    verify_guards(guards)
    if hashlib.sha256(checkpoint_path.read_bytes()).hexdigest() != checkpoint_sha:
        raise RuntimeError('Finite checkpoint changed.')
    write_packet(out, result)
    print(verdict, 'seconds', result['elapsed_seconds'], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=HERE)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.with_suffix('.finite.zlib').exists() or args.out.suffix != '.zlib':
        parser.error('Use a fresh .zlib output; old archives are never overwritten.')
    ctx.prec = 1024
    consume(args.input_dir.resolve(), args.out)


if __name__ == '__main__':
    main()
