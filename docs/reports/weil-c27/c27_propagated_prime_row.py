"""Certify one canonical propagated weighted row of the actual c27 prime operator.

Reuse positive weighted Schur and quadratic Bernstein bounds. Only the
ordered two-step clipping geometry is new; accepted producers are not run.
The supplied source norm and excluded coefficient endpoint are read as data.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from math import isqrt
from pathlib import Path
import time
import zlib

from flint import arb, ctx, fmpq
from joint_gamma_boundary import write_packet
from structured_inverse_residual_interval import envelope

HERE = Path(__file__).resolve().parent
BITS = 640
SOURCE_ID = '1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e'
PINS = {
    'high-floor.json': '3e85480a54d57e2e842964d9be764c755aedfeefd3ff6feec3f851dc0acf8120',
    'full-gamma-surplus-budget.json': '25eedbaa0dc4fe6612b58cd34f4fbf544719e5d373c2fbccd4b2d11a8b64597f',
    'full-gamma-stationarity-family-norm.json.zlib': '169e858deb9f06ca677098da032849450074f794627eb0174711717abc5ca1b8',
}
HELPERS = {
    'joint_gamma_boundary.py': '10021a4faf837683e8430a7124241f3145e1bcfdce42aca3f817cdd4a491ea07',
    'structured_inverse_residual_interval.py': 'eae66659370f74c8a4c7adbf8df59755c8aa25d532a9652a3c16e7f2da10567c',
}
PRIME_POWERS = {2:(2,1), 3:(3,1), 4:(2,2), 5:(5,1), 7:(7,1),
    8:(2,3), 9:(3,2), 11:(11,1), 13:(13,1), 16:(2,4), 17:(17,1),
    19:(19,1), 23:(23,1), 25:(5,2)}


def ball(q):
    q = F(q)
    return arb(fmpq(q.numerator, q.denominator))


def rat(q):
    q = F(q)
    return [str(q.numerator), str(q.denominator)]


def switch_labels(shifts):
    """e=exp(L-x): an ordered path needs a<e<27a and ab<e<27ab."""
    factors = set(shifts)
    factors.update(a*b for a in shifts for b in shifts)
    events = factors | {27*a for a in factors}
    return sorted((a for a in events if 1 < a and a*a < 27), reverse=True)


def cell_sample(left, right):
    """left=None denotes sqrt(27); sample only selects the exact active set."""
    top = F(isqrt(27*10**36), 10**18) if left is None else left
    if not top > right:
        raise ArithmeticError('No certified rational interior sample.')
    return (top+right)/2


def bernstein_term(shift, left, right):
    """Degree2 Bernstein coefficients on the exact, possibly uncertain cell."""
    a, b = left+shift, right+shift
    c = ball(F(29,50))
    return (1+c*a*a, 1+c*a*b, 1+c*b*b)


def ratio_upper(numerator, denominator):
    if (len(numerator) != 3 or len(denominator) != 3
            or any(not x.is_finite() for x in (*numerator, *denominator))
            or any(not d > 0 for d in denominator)):
        raise ArithmeticError('Missing positive Bernstein denominator.')
    # With nonnegative Bernstein basis, N/D is a convex combination of N_i/D_i.
    return max((n/d).upper() for n,d in zip(numerator, denominator))


def propagated_row(cache):
    logs = {}
    for p, bounds in cache['logarithm_intervals'].items():
        lo, hi = map(F, bounds)
        if not 0 < lo <= hi:
            raise ValueError('Invalid accepted logarithm enclosure.')
        logs[int(p)] = ball((lo+hi)/2)+arb(0, ball((hi-lo)/2))
    if set(logs) != {p for p,_ in PRIME_POWERS.values()}:
        raise ValueError('Incomplete prime logarithms.')
    width = 3*logs[3]/2
    if not width > 0:
        raise ArithmeticError('Missing positive half-width.')

    @lru_cache(None)
    def normalized_log(q):
        numerator, denominator = q.numerator, q.denominator
        value = arb(0)
        for p, log in logs.items():
            exponent = 0
            while numerator % p == 0:
                numerator //= p
                exponent += 1
            while denominator % p == 0:
                denominator //= p
                exponent -= 1
            value += exponent*log
        if numerator != 1 or denominator != 1:
            raise ValueError('Unrepresented rational switch factor.')
        return value/width

    shifts = {}
    for n, (p,k) in PRIME_POWERS.items():
        if p**k != n:
            raise ValueError('Invalid prime-power weight metadata.')
        weight = logs[p]/arb(n).sqrt()
        if not weight > 0:
            raise ArithmeticError('Nonpositive prime-power weight.')
        for a in (F(n), F(1,n)):
            shifts[a] = (normalized_log(a), weight)
    labels = [None]+switch_labels(shifts)+[F(1)]
    rows = []
    for left_label, right_label in zip(labels, labels[1:]):
        sample = cell_sample(left_label, right_label)
        left = arb(0) if left_label is None else 1-normalized_log(left_label)
        right = 1-normalized_log(right_label)
        if not right-left > 0:
            raise ArithmeticError('Unresolved exact support-cell order.')
        first = [a for a in shifts if a < sample < 27*a]
        if F(1,2) not in first:
            raise ArithmeticError('Missing everywhere-active negative log2 shift.')
        denominator, numerator = [arb(0) for _ in range(3)], [arb(0) for _ in range(3)]
        paths = 0
        for a in first:
            shift, weight = shifts[a]
            for k,v in enumerate(bernstein_term(shift,left,right)):
                denominator[k] += weight*v
            for b, (_,weight_b) in shifts.items():
                ab = a*b
                if ab < sample < 27*ab:
                    paths += 1
                    for k,v in enumerate(bernstein_term(normalized_log(ab),left,right)):
                        numerator[k] += weight*weight_b*v
        upper = ratio_upper(numerator, denominator)
        exported = {key:[envelope(v,BITS) for v in values] for key,values in
            dict(denominator_bernstein=denominator,numerator_bernstein=numerator).items()}
        rows.append(dict(left_switch=None if left_label is None else rat(left_label),
            right_switch=rat(right_label), normalized_left=envelope(left,BITS),
            normalized_right=envelope(right,BITS), active_first=[rat(a) for a in first],
            active_ordered_paths=paths, row_upper=envelope(arb(upper),BITS), **exported))
    if not rows:
        raise ArithmeticError('Empty clipped-path cover.')
    upper = max(F(int(row['row_upper'][1]),1 << BITS) for row in rows)
    return upper, rows


def consume(out):
    started = time.perf_counter()
    packets = {}
    for name, digest in dict(PINS, **HELPERS).items():
        raw = (HERE/name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('Different accepted input or export helper: '+name)
        if name in PINS:
            packets[name] = json.loads(zlib.decompress(raw) if name.endswith('.zlib') else raw)
    cache, budget, endpoint = (packets[n] for n in PINS)
    if ((cache['cutoff'],cache['half_width'],cache['weight'],cache['row_upper']) !=
            (27,'3log(3)/2','1+(29/50)(u/a)^2','59/10')
            or {int(n):tuple(pk) for n,pk in cache['prime_powers'].items()} != PRIME_POWERS
            or budget['schema'] != 'ACTUAL_C27_CACHED_SURPLUS_BUDGET_V1'
            or endpoint['schema'] != 'ACTUAL_C27_STATIONARITY_COHERENT_FAMILY_V1'
            or any((p['source_integer_sha256'],p['source_bits'],p['source_row']) !=
                   (SOURCE_ID,384,0) for p in (budget,endpoint))):
        raise ValueError('Different actual row, SAME source or coefficient endpoint.')
    get = lambda p,k: F(*map(int,p['exact_rationals'][k]))
    norm = get(budget,'full_unweighted_norm_squared_upper')
    hlow, hupp = get(budget,'h_lower'), get(budget,'h_upper')
    old_coefficient = get(endpoint,'t_plus')
    if not norm > 0 or not 0 < hlow <= hupp or not old_coefficient > 0:
        raise ArithmeticError('Missing complete norm, positive h or coefficient cap.')
    upper, rows = propagated_row(cache)
    gap = F(59,10)-upper
    values = dict(propagated_row_upper=upper, actual_surplus_floor_lower=gap,
        complete_full_source_norm_squared_upper=norm, h_lower=hlow, h_upper=hupp,
        inherited_excluded_rank_one_coefficient_endpoint=old_coefficient)
    if gap > 0:
        rho, tau = gap/norm, gap*hlow/norm
        values.update(directional_budget_upper=norm/gap,
            certified_rank_one_coefficient=rho, normalized_tau_lower=tau,
            shrink_parameter_lower=tau/(1+tau),
            coefficient_gain_over_excluded_endpoint=rho-old_coefficient)
    verdict = ('actual_slack_and_new_eligible_coefficient' if gap > 0 and
        values['certified_rank_one_coefficient'] > old_coefficient else
        'improved_actual_slack_but_insufficient_coefficient' if gap > F(1,632) else
        'stopped_canonical_propagated_row_insufficient')
    for name, digest in dict(PINS, **HELPERS).items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise RuntimeError('Accepted input changed during propagation: '+name)
    result = dict(schema='ACTUAL_C27_CANONICAL_PROPAGATED_PRIME_ROW_V1',
        packet_sha256=PINS, inherited_dependency_sha256=HELPERS,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        grid_bits=BITS, working_bits=ctx.prec, source_integer_sha256=SOURCE_ID,
        source_bits=384, source_row=0, canonical_weight='q1=C_p*(1+(29/50)*(x/L)^2)',
        source_support=[2049,32768], complete_far_retained=True,
        exact_rationals={k:rat(v) for k,v in values.items()}, cell_results=rows,
        support_cells=len(rows), ordered_path_cell_count=sum(r['active_ordered_paths'] for r in rows),
        verdict=verdict, old_scientific_producers_or_fixtures_replayed=0,
        new_source_actions_Grams_inverses_solves_or_far_queries=0,
        propagation_applications=1, propagation_order_weight_parameter_precision_or_partition_scan=False,
        actual_retained_or_cofinal_positivity_claim=False, actual_rh_proof=False,
        elapsed_seconds=time.perf_counter()-started, conditions=budget['conditions']+[
            'ALL actual common-form/domain and original source/normalization/cache identifications remain premises.',
            'C_p is the actual selfadjoint clipped prime-power sum on I=(-L,L). q1=C_p*q0 is bounded and strictly positive; the negative log2 shift is active on0<=x<=L. Reflection covers the other half.',
            'Each exact support cell includes every ordered path with both intermediate and final points in I. Exact rational switch ordering supplies active sets; isolated endpoints have zero measure.',
            'Positive denominator Bernstein coefficients give C_p*q1/q1 <= max_i(N_i/D_i) on each cell. Outward coefficient and endpoint intervals cover the exact quantities.',
            'Weighted Schur gives C_p<=R1 I and Hext-T>=(59/10-R1)I. Infimum over the SAME legal finite E shifts gives Pi>=gap I on O; no band invariance is required.',
            'The saved complete full norm includes every far term. uu*<=N I<=(N/gap)Pi implies M+(gap/N)uu*<=M+Pi=A. tau=h*gap/N,while the excluded tplus is a rank-one coefficient; these scales are not identified.',
            'An eligible coefficient is not a sign of the actual retained Weil form,cofinal positivity,Robin or RH. Standard Schur/Bernstein methods are reused; no generic novelty or Lean claim.',
            'One canonical propagated row only. No iteration,source,weight,parameter,precision or partition retry. The previously excluded lower-comparison route remains stopped.'])
    write_packet(out,result)
    print(verdict,'cells',len(rows),'seconds',result['elapsed_seconds'],flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.suffix != '.zlib':
        parser.error('Use one fresh .zlib output; no canonical propagation retry.')
    ctx.prec = 1024
    consume(args.out)


if __name__ == '__main__':
    main()
