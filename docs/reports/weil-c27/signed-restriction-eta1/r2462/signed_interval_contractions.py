"""Candidate fixed W1 contraction supplier; no input acquisition or CLI.

Public replay fixes the original authenticated1263 modes,133 RHS, eta1,
1024-bit arithmetic, accepted plan and accepted reference transcript.
The parameterized private loop is for small exact dense validation only.
No large target run has been performed to prepare this candidate.
"""
from fractions import Fraction
import hashlib
from pathlib import Path

CORE_SHA256 = 'b278f8a358845ca436166b5595f8043173c7e09678d805d9d2c7e54d9d83d1c7'
BASE_SHA256 = '8cdfce25b0342e0b0f5b208f186aef1633e69209383a178e42311b76c689256e'
PLAN_SHA256 = 'dc010d3b6ff069d67bf168786a275806233d49adc7be74710329cf945019b8b3'
REFERENCE_SHA256 = '2a0f2e748634c355d24624e5dbc4b1f87da3afca4d5651894c235cb7665544ac'
DIMENSION, FEATURES, RHS_WIDTH, ARITHMETIC_BITS = 1263, 131, 133, 1024


def _zero(core, bits):
    return core.Interval(0, 0, bits)


def _transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def _matmul(core, left, right):
    core.require(bool(left) and bool(right) and bool(right[0]), 'shape_mismatch', 'matrix product')
    inner, cols = len(right), len(right[0])
    core.require(all(len(row) == inner for row in left)
                 and all(len(row) == cols for row in right), 'shape_mismatch', 'matrix product')
    bits = left[0][0].bits
    result = []
    for row in left:
        out = []
        for j in range(cols):
            value = _zero(core, bits)
            for k in range(inner):
                value = value + row[k] * right[k][j]
            out.append(value)
        result.append(out)
    return result


def _symmetrize(core, matrix):
    size = len(matrix)
    core.require(all(len(row) == size for row in matrix), 'shape_mismatch', 'symmetric square')
    for i in range(size):
        for j in range(i + 1, size):
            shared = matrix[i][j].intersect(matrix[j][i])
            matrix[i][j] = matrix[j][i] = shared
    return matrix


def _ordinary_contractions(core, labels, features, actions):
    """Unshifted F^T X and X^T X; no inverse is computed."""
    count = len(features[labels[0]])
    bits = features[labels[0]][0].bits
    mixed = [[_zero(core, bits) for _ in range(2)] for _ in range(count)]
    euclidean = [[_zero(core, bits) for _ in range(2)] for _ in range(2)]
    for n in labels:
        for i in range(count):
            for j in range(2):
                mixed[i][j] = mixed[i][j] + features[n][i] * actions[n][j]
        euclidean[0][0] = euclidean[0][0] + actions[n][0].square()
        euclidean[0][1] = euclidean[0][1] + actions[n][0] * actions[n][1]
        euclidean[1][1] = euclidean[1][1] + actions[n][1].square()
    euclidean[1][0] = euclidean[0][1]
    return mixed, euclidean


def _replay_state_with_rhs(core, base, labels, diagonal, generators, eta_sigma,
                           blocks, dimension, rhs, width, reference=None):
    """The accepted matrix loop plus forward RHS and inverse-Gram updates.

    dimensions/width may vary only for private synthetic tests; the public
    authenticated caller always supplies1263 and133 and the fixed reference.
    """
    bits = eta_sigma.bits
    core.require(type(width) is int and width > 0 and len(labels) == dimension
                 and len(set(labels)) == dimension and set(rhs) == set(labels)
                 and all(len(rhs[n]) == width and all(x.bits == bits for x in rhs[n])
                         for n in labels), 'shape_mismatch', 'RHS state')
    if reference is not None:
        core.require(len(reference) == len(blocks), 'reference_mismatch', 'reference step count')
    active = list(labels)
    positive = negative = 0
    transcript = []
    gram = [[_zero(core, bits) for _ in range(width)] for _ in range(width)]
    for step, pivot_labels in enumerate(blocks):
        core.require(all(n in active for n in pivot_labels),
                     'invalid_plan', 'a pivot label was already removed')
        rest = [n for n in active if n not in pivot_labels]
        pivot_rows = [list(generators[n]) for n in pivot_labels]
        if len(pivot_labels) == 1:
            k = pivot_labels[0]
            p = diagonal[k]
            inverse = [[base._reciprocal(core, p, f'step {step}, label {k}')]]
            signs = (1, 0) if p.lo > 0 else (0, 1)
            pivot_matrix = [[p.export()]]
            determinant = p
        else:
            k, l = pivot_labels
            a, d = diagonal[k], diagonal[l]
            b = base._entry(core, k, l, generators, eta_sigma)
            determinant = a * d - b.square()
            invdet = base._reciprocal(core, determinant,
                                      f'step {step}, labels {k},{l}, determinant')
            inverse = [[d * invdet, -(b * invdet)],
                       [-(b * invdet), a * invdet]]
            if determinant.hi < 0:
                signs = (1, 1)
            else:
                trace = a + d
                core.require(trace.lo > 0 or trace.hi < 0,
                             'uncertain_pivot_sign',
                             f'step {step}, labels {k},{l}: positive determinant, uncertain trace')
                signs = (2, 0) if trace.lo > 0 else (0, 2)
            pivot_matrix = [[a.export(), b.export()], [b.export(), d.export()]]
        size = len(pivot_labels)
        # New work: snapshot pre-update RHS and add Q_K^T A^-1 Q_K.
        pivot_rhs = [list(rhs[n]) for n in pivot_labels]
        for a in range(width):
            for b in range(a, width):
                contribution = _zero(core, bits)
                for u in range(size):
                    for v in range(size):
                        product = (pivot_rhs[u][a].square() if a == b and u == v
                                   else pivot_rhs[u][a] * pivot_rhs[v][b])
                        contribution = contribution + product * inverse[u][v]
                gram[a][b] = gram[a][b] + contribution
        for i in rest:
            # Original matrix work: all entries use the old current generators.
            c = [base._entry(core, i, k, generators, eta_sigma) for k in pivot_labels]
            multipliers = []
            for j in range(size):
                value = core.Interval(0, 0, bits)
                for t in range(size):
                    value = value + c[t] * inverse[t][j]
                multipliers.append(value)
            # New work: use exactly these same pre-update multipliers.
            for column in range(width):
                value = rhs[i][column]
                for j in range(size):
                    value = value - multipliers[j] * pivot_rhs[j][column]
                rhs[i][column] = value
            debit = core.Interval(0, 0, bits)
            for j in range(size):
                debit = debit + multipliers[j] * c[j]
            diagonal[i] = diagonal[i] - debit
            old_row = list(generators[i])
            for column in range(4):
                value = old_row[column]
                for j in range(size):
                    value = value - multipliers[j] * pivot_rows[j][column]
                generators[i][column] = value
        positive += signs[0]
        negative += signs[1]
        record = {'step': step, 'original_labels': list(pivot_labels),
                  'pivot_interval_matrix': pivot_matrix,
                  'determinant_interval': determinant.export(),
                  'inertia': [signs[0], signs[1], 0]}
        if reference is not None:
            core.require(record == reference[step], 'reference_mismatch', f'pivot step {step}')
        transcript.append(record)
        active = rest
    core.require(not active and positive + negative == dimension,
                 'incomplete_coverage', 'no final count or Gram from a prefix')
    for a in range(width):
        for b in range(a + 1, width):
            gram[b][a] = gram[a][b]
    return positive, negative, transcript, gram


def _recover_border(core, base, gram, mixed, euclidean, saved_metric, kernel, epsilon, far):
    """Small interval contractions only. No kernel/block inverse or sign claim."""
    features = len(kernel)
    bits = epsilon.bits
    core.require(0 < far < features and len(gram) == features + 2
                 and all(len(row) == features + 2 for row in gram),
                 'shape_mismatch', 'small contraction blocks')
    h = [row[:features] for row in gram[:features]]
    hx = [row[features:] for row in gram[:features]]
    xx = [row[features:] for row in gram[features:]]
    ell = _matmul(core, kernel, mixed)  # L=Kt F^T X; not the combined kernel.
    hl = _matmul(core, h, ell)
    j = [[mixed[a][b] + epsilon * hx[a][b] - hl[a][b]
          for b in range(2)] for a in range(features)]
    atl = _matmul(core, _transpose(mixed), ell)
    hxtl = _matmul(core, _transpose(hx), ell)
    lthx = _matmul(core, _transpose(ell), hx)
    lthl = _matmul(core, _transpose(ell), hl)
    epsilon2 = epsilon.square()
    paid = [[saved_metric[a][b] + epsilon * euclidean[a][b] - atl[a][b]
             + epsilon2 * xx[a][b] - epsilon * (hxtl[a][b] + lthx[a][b])
             + lthl[a][b] for b in range(2)] for a in range(2)]
    paid = _symmetrize(core, paid)
    combined = []
    for a in range(features):
        row = []
        for b in range(features):
            if a < far and b < far:
                value = base._scale_rational(core, kernel[a][b], Fraction(65, 33))
            elif a >= far and b >= far:
                value = kernel[a][b]
            else:
                core.require(kernel[a][b].is_zero(), 'kernel_structure_mismatch', 'cross block')
                value = _zero(core, bits)
            row.append(value)
        combined.append(row)
    khk = _matmul(core, _matmul(core, combined, h), combined)
    kj = _matmul(core, combined, j)
    border = [[_zero(core, bits) for _ in range(features + 2)] for _ in range(features + 2)]
    for a in range(features):
        for b in range(features):
            border[a][b] = combined[a][b] - khk[a][b]
        for b in range(2):
            border[a][features + b] = border[features + b][a] = -kj[a][b]
    for a in range(2):
        for b in range(2):
            border[features + a][features + b] = -paid[a][b]
    border = _symmetrize(core, border)
    return j, paid, combined, border


def _export_matrix(matrix):
    return [[entry.export() for entry in row] for row in matrix]


def _export_upper(matrix):
    return [[matrix[i][j].export() for j in range(i, len(matrix))]
            for i in range(len(matrix))]


def replay(core, base, target_raw, metric_raw, plan_raw, reference_raw):
    """One fixed source1263/133 research pass; caller owns resources and IO."""
    core.require(hashlib.sha256(Path(core.__file__).read_bytes()).hexdigest() == CORE_SHA256,
                 'core_source_mismatch', 'pinned interval core required')
    core.require(hashlib.sha256(Path(base.__file__).read_bytes()).hexdigest() == BASE_SHA256,
                 'base_source_mismatch', 'pinned executed r2429 source required')
    for raw, expected, where in ((plan_raw, PLAN_SHA256, 'plan'),
                                  (reference_raw, REFERENCE_SHA256, 'accepted result')):
        core.require(type(raw) is bytes and hashlib.sha256(raw).hexdigest() == expected,
                     'binding_mismatch', where + ' fixed immutable bytes')
    plan, reference = core.parse_json(plan_raw), core.parse_json(reference_raw)
    for key, expected in {
        'schema': base.PLAN_SCHEMA, 'purpose': 'inertia_research_only',
        'target_sha256': core.TARGET_SHA256, 'metric_decoded_sha256': core.METRIC_SHA256,
        'source_integer_sha256': core.SOURCE_ID,
        'source_structure_proof_sha256': base.SOURCE_PROOF_SHA256,
        'arithmetic_bits': ARITHMETIC_BITS, 'eta': ['1', '1'],
    }.items():
        core.equal_field(plan, key, expected, 'plan')
    eta = base._rational(core, plan['eta'], 'plan.eta')
    bits = core.integer(plan['arithmetic_bits'], 'plan.arithmetic_bits')
    blocks = plan['blocks']
    core.require(type(blocks) is list and len(blocks) == DIMENSION
                 and all(type(b) is list and len(b) == 1 and type(b[0]) is int for b in blocks),
                 'invalid_plan', 'fixed1263 singleton plan')
    flat = [b[0] for b in blocks]
    core.require(len(set(flat)) == DIMENSION, 'incomplete_coverage', 'distinct fixed original labels')
    labels, remaining, ell, full_generators, sigma, terms = base._decode_bound_inputs(
        core, target_raw, metric_raw)
    core.require(flat == labels, 'partition_mismatch', 'fixed saved order required')
    for key, expected in {'schema': base.RESULT_SCHEMA,
                          'status': 'complete_inertia_conditional_on_inherited_source_contract',
                          'eta': ['1', '1'], 'arithmetic_bits': bits,
                          'inertia': [DIMENSION, 0, 0], 'original_good_modes': labels,
                          'original_remaining_modes': remaining}.items():
        core.equal_field(reference, key, expected, 'accepted_reference')
    core.equal_field(reference['input_bindings'], 'plan_sha256', PLAN_SHA256, 'accepted_reference')
    metric = core.parse_json(metric_raw)  # The exact bytes were authenticated by the base decoder.
    raw_x = core.decode_matrix(metric['inverse_actions'], DIMENSION, 2, 'metric.inverse_actions')
    saved = core.decode_matrix(metric['joint_inverse_metric'], 2, 2, 'metric.joint_inverse_metric')
    saved = _symmetrize(core, [[base._lift(core, x, bits) for x in row] for row in saved])
    core.require(saved[0][0].lo > 0 and (saved[0][0] * saved[1][1] - saved[0][1].square()).lo > 0,
                 'uncertain_paid_rank', 'saved same-source joint metric must certify rank two')
    epsilon = core.Interval.rational(
        Fraction(27852372618143874844, 1 << 512) + Fraction(1, 1 << 542), bits=bits)
    eta_sigma = base._scale_rational(core, base._lift(core, sigma, bits), eta)
    diagonal, generators, features, actions = {}, {}, {}, {}
    for position, n in enumerate(labels):
        c_diagonal = ell[n] + core.quadratic(full_generators[n][4:4 + core.RANK], terms)
        diagonal[n] = base._scale_rational(core, base._lift(core, c_diagonal, bits), eta) - epsilon
        generators[n] = [base._lift(core, x, bits) for x in full_generators[n][:4]]
        features[n] = [base._lift(core, x, bits) for x in full_generators[n][4:4 + FEATURES]]
        actions[n] = [base._lift(core, x, bits) for x in raw_x[position]]
    kernel = [[_zero(core, bits) for _ in range(FEATURES)] for _ in range(FEATURES)]
    for i, j, coefficient in terms:
        kernel[i][j] = kernel[j][i] = base._lift(core, coefficient, bits)
    del ell, full_generators, terms, raw_x, metric
    mixed, euclidean = _ordinary_contractions(core, labels, features, actions)
    rhs = {n: list(features[n]) + list(actions[n]) for n in labels}
    del features, actions
    positive, negative, transcript, gram = _replay_state_with_rhs(
        core, base, labels, diagonal, generators, eta_sigma, blocks, DIMENSION,
        rhs, RHS_WIDTH, reference['pivot_transcript'])
    core.require((positive, negative) == (DIMENSION, 0), 'reference_mismatch', 'positive fixed base')
    j, paid, combined, border = _recover_border(
        core, base, gram, mixed, euclidean, saved, kernel, epsilon, core.FAR)
    return {
        'schema': 'R2462_FIXED_W1_FORWARD_GRAM_V1',
        'status': 'complete_contractions_conditional_on_inherited_source_and_saved_action_contract',
        'eta': ['1', '1'], 'arithmetic_bits': bits,
        'dimension': DIMENSION, 'rhs_width': RHS_WIDTH, 'feature_count': FEATURES,
        'rhs_order': 'first131 original F columns, then the two original metric X columns',
        'base_inertia': [positive, negative, 0], 'pivot_transcript': transcript,
        'all_pivots_equal_accepted_reference': True,
        'original_good_modes': labels, 'original_remaining_modes': remaining,
        'symmetric_encoding': 'upper row i stores columns i through dimension-1',
        'gram_F_X_upper': _export_upper(gram),
        'F_transpose_X': _export_matrix(mixed), 'X_transpose_X_upper': _export_upper(euclidean),
        'saved_metric_G_upper': _export_upper(saved), 'original_Kt_upper': _export_upper(kernel),
        'combined_kernel_C1_upper': _export_upper(combined),
        'F_transpose_Winv_Btranspose': _export_matrix(j), 'B_Winv_Btranspose_upper': _export_upper(paid),
        'inverse_free_border133_upper': _export_upper(border),
        'input_bindings': dict(reference['input_bindings'],
                               reference_result_sha256=REFERENCE_SHA256,
                               matrix_replay_source_sha256=BASE_SHA256,
                               contraction_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),
        'padding_account': dict(reference['padding_account']),
        'premises': list(reference['premises']) + [
            'Saved same-source PX=B^T and G=BX are inherited action identities, not proved by a zero-containing residual.',
            'Only Gram contractions are computed; no full inverse or inverse-action vectors are materialized.',
            'The inverse-free133 matrix has not received a signed inertia certificate.'],
        'border_inertia': None, 'restricted_domination_claim': False,
        'eta_budget_membership_claim': False, 'tau_or_z_upper_claim': False,
        'MAIN_or_RH_claim': False, 'compiled_or_official_acceptance': False,
    }
