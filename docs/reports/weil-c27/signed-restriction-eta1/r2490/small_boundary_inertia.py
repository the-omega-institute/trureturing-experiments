"""Fixed paid2/boundary19 reduction and prescribed small signed inertia.

No IO/acquisition/CLI, pivot search, target replay, or default plan. The caller
provides immutable result bytes and their independently nominated SHA256.
The private dense kernel is parameterized only for synthetic exact checks.
All conclusions remain conditional on the original common-source premises.
"""
import hashlib
from pathlib import Path


CORE_SHA256 = 'b278f8a358845ca436166b5595f8043173c7e09678d805d9d2c7e54d9d83d1c7'
BASE_SHA256 = '8cdfce25b0342e0b0f5b208f186aef1633e69209383a178e42311b76c689256e'
CONTRACTION_SHA256 = '2f54ee32dcbcbc934f31c8b9f10772322589ae79c8164071f432bd0321785b37'
ORIGINAL_PLAN_SHA256 = 'dc010d3b6ff069d67bf168786a275806233d49adc7be74710329cf945019b8b3'
REFERENCE_SHA256 = '2a0f2e748634c355d24624e5dbc4b1f87da3afca4d5651894c235cb7665544ac'
FAR, BOUNDARY, PAID, BITS = 112, 19, 2, 1024
FEATURES, BORDER = FAR + BOUNDARY, FAR + BOUNDARY + PAID
CONTRACTION_SCHEMA = 'R2462_FIXED_W1_FORWARD_GRAM_V1'
REDUCTION_SCHEMA = 'R2472_FIXED_PAID2_BOUNDARY19_REDUCTION_V1'
PLAN_SCHEMA = 'R2472_PRESCRIBED_REDUCED112_PLAN_V1'
RESULT_SCHEMA = 'R2472_REDUCED112_SIGNED_INERTIA_V1'
ENCODING = 'upper row i stores columns i through dimension-1'
PADDING = {
    'C_initialization': 'ell plus the complete shared131-feature quadratic',
    'embedded_charge': 'a_orig + (1/32)*a_reuse; neither recovered separately nor altered',
    'additional_shift': 'eta*C - epsilon*I; complete epsilon charged once',
    'epsilon_exact': '27852372618143874844/2^512 + 1/2^542',
}


def _zero(core, bits):
    return core.Interval(0, 0, bits)


def _authenticate(core, raw, expected_sha256, where):
    core.require(type(expected_sha256) is str and len(expected_sha256) == 64
                 and all(c in '0123456789abcdef' for c in expected_sha256),
                 'invalid_expected_hash', where + ': caller-nominated SHA256 required')
    core.require(type(raw) is bytes
                 and hashlib.sha256(raw).hexdigest() == expected_sha256,
                 'binding_mismatch', where + ': exact immutable bytes')
    return core.parse_json(raw)


def _sources(core, base):
    core.require(hashlib.sha256(Path(core.__file__).read_bytes()).hexdigest() == CORE_SHA256,
                 'core_source_mismatch', 'pinned outward interval core')
    core.require(hashlib.sha256(Path(base.__file__).read_bytes()).hexdigest() == BASE_SHA256,
                 'base_source_mismatch', 'pinned signed reciprocal source')


def _bindings(core, base, value):
    core.require(type(value) is dict, 'binding_mismatch', 'input_bindings object')
    for key, expected in {
        'target_git_blob': core.TARGET_BLOB,
        'target_sha256': core.TARGET_SHA256,
        'metric_decoded_sha256': core.METRIC_SHA256,
        'source_integer_sha256': core.SOURCE_ID,
        'basis': core.BASIS,
        'core_source_sha256': CORE_SHA256,
        'source_structure_proof_sha256': base.SOURCE_PROOF_SHA256,
        'plan_sha256': ORIGINAL_PLAN_SHA256,
        'reference_result_sha256': REFERENCE_SHA256,
        'matrix_replay_source_sha256': BASE_SHA256,
        'contraction_source_sha256': CONTRACTION_SHA256,
    }.items():
        core.equal_field(value, key, expected, 'input_bindings')


def _decode_upper(core, raw, size, bits, where):
    core.require(type(raw) is list and len(raw) == size, 'shape_mismatch', where)
    matrix = [[None] * size for _ in range(size)]
    for i, row in enumerate(raw):
        core.require(type(row) is list and len(row) == size - i,
                     'shape_mismatch', f'{where}[{i}]')
        for offset, pair in enumerate(row):
            j = i + offset
            value = core.Interval.decode(pair, f'{where}[{i}][{offset}]', bits=bits)
            matrix[i][j] = matrix[j][i] = value
    return matrix


def _export_upper(matrix):
    return [[matrix[i][j].export() for j in range(i, len(matrix))]
            for i in range(len(matrix))]


def _check_blocks(core, labels, blocks, complete=False):
    core.require(type(blocks) is list and bool(blocks), 'invalid_plan', 'explicit blocks required')
    flat = []
    for block in blocks:
        core.require(type(block) is list and len(block) in (1, 2)
                     and all(type(i) is int for i in block),
                     'invalid_plan', 'original index blocks have size one or two')
        flat.extend(block)
    core.require(len(flat) == len(set(flat)) and set(flat) <= set(labels),
                 'invalid_plan', 'distinct covered original indices required')
    if complete:
        core.require(len(flat) == len(labels) and set(flat) == set(labels),
                     'incomplete_coverage', 'no complete inertia from a prefix')


def _eliminate(core, base, matrix, labels, blocks, required_inertia=None):
    """Enclose the Schur complement for prescribed original-label blocks.

    Returns the remaining matrix/labels and *removed* inertia, never an inertia
    of an unfinished suffix. Input intervals are copied and not modified.
    """
    size = len(labels)
    core.require(size > 0 and len(set(labels)) == size
                 and all(type(i) is int for i in labels)
                 and len(matrix) == size and all(len(row) == size for row in matrix),
                 'shape_mismatch', 'dense symmetric state')
    bits = matrix[0][0].bits
    core.require(all(type(x) is core.Interval and x.bits == bits for row in matrix for x in row)
                 and all(matrix[i][j] == matrix[j][i] for i in range(size) for j in range(size)),
                 'invalid_symmetric_state', 'one common grid and mirrored intervals')
    _check_blocks(core, labels, blocks)
    if required_inertia is not None:
        core.require(len(required_inertia) == len(blocks), 'invalid_plan', 'required sign count')
    active, current = list(labels), [list(row) for row in matrix]
    positive = negative = 0
    transcript = []
    for step, block in enumerate(blocks):
        positions = {label: i for i, label in enumerate(active)}
        pivot_positions = [positions[i] for i in block]
        rest = [i for i in active if i not in block]
        rest_positions = [positions[i] for i in rest]
        if len(block) == 1:
            p = current[pivot_positions[0]][pivot_positions[0]]
            inverse = [[base._reciprocal(core, p, f'dense step {step}, index {block[0]}')]]
            signs = [1, 0, 0] if p.lo > 0 else [0, 1, 0]
            pivot = [[p]]
            determinant = p
        else:
            k, l = pivot_positions
            a, b, d = current[k][k], current[k][l], current[l][l]
            determinant = a * d - b.square()
            invdet = base._reciprocal(core, determinant, f'dense step {step}, determinant')
            inverse = [[d * invdet, -(b * invdet)], [-(b * invdet), a * invdet]]
            if determinant.hi < 0:
                signs = [1, 1, 0]
            else:
                trace = a + d
                core.require(trace.lo > 0 or trace.hi < 0, 'uncertain_pivot_sign',
                             f'dense step {step}: positive determinant, uncertain trace')
                signs = [2, 0, 0] if trace.lo > 0 else [0, 2, 0]
            pivot = [[a, b], [b, d]]
        if required_inertia is not None:
            core.require(signs == required_inertia[step], 'wrong_pivot_sign',
                         f'dense step {step}, original indices {block}')
        width = len(block)
        cross = [[current[i][k] for k in pivot_positions] for i in rest_positions]
        multipliers = []
        for row in cross:
            out = []
            for j in range(width):
                value = _zero(core, bits)
                for k in range(width):
                    value = value + row[k] * inverse[k][j]
                out.append(value)
            multipliers.append(out)
        next_matrix = [[None] * len(rest) for _ in rest]
        for i, old_i in enumerate(rest_positions):
            for j in range(i, len(rest)):
                value = current[old_i][rest_positions[j]]
                for k in range(width):
                    value = value - multipliers[i][k] * cross[j][k]
                next_matrix[i][j] = next_matrix[j][i] = value
        transcript.append({
            'step': step, 'original_indices': list(block),
            'pivot_interval_matrix': [[x.export() for x in row] for row in pivot],
            'determinant_interval': determinant.export(), 'inertia': signs,
        })
        positive += signs[0]
        negative += signs[1]
        active, current = rest, next_matrix
    core.require(positive + negative + len(active) == size,
                 'incomplete_coverage', 'each pivot counted exactly once')
    return current, active, [positive, negative, 0], transcript


def _reduce_fixed(core, base, matrix, far, boundary):
    """Private size parameters support small same-block synthetic fixtures."""
    core.require(type(far) is int and far > 0 and type(boundary) is int and boundary > 0
                 and len(matrix) == far + boundary + PAID,
                 'shape_mismatch', 'fixed far/boundary/paid ordering')
    features = far + boundary
    paid_blocks = [[features], [features + 1]]
    boundary_blocks = [[i] for i in range(far, features)]
    after_paid, feature_labels, paid_inertia, paid_transcript = _eliminate(
        core, base, matrix, list(range(features + PAID)), paid_blocks, [[0, 1, 0]] * PAID)
    reduced, retained, boundary_inertia, boundary_transcript = _eliminate(
        core, base, after_paid, feature_labels, boundary_blocks, [[0, 1, 0]] * boundary)
    core.require(retained == list(range(far)) and paid_inertia == [0, PAID, 0]
                 and boundary_inertia == [0, boundary, 0], 'incomplete_coverage', 'fixed reduction')
    return after_paid, reduced, paid_transcript, boundary_transcript


def reduce_paid_boundary(core, base, contraction_raw, expected_contraction_sha256):
    """Only the fixed original133 ->131 ->112 reduction; no final sign count."""
    _sources(core, base)
    packet = _authenticate(core, contraction_raw, expected_contraction_sha256, 'contraction result')
    for key, expected in {
        'schema': CONTRACTION_SCHEMA,
        'status': 'complete_contractions_conditional_on_inherited_source_and_saved_action_contract',
        'dimension': 1263, 'rhs_width': BORDER, 'feature_count': FEATURES,
        'rhs_order': 'first131 original F columns, then the two original metric X columns',
        'arithmetic_bits': BITS, 'eta': ['1', '1'], 'base_inertia': [1263, 0, 0],
        'all_pivots_equal_accepted_reference': True, 'symmetric_encoding': ENCODING,
        'padding_account': PADDING,
    }.items():
        core.equal_field(packet, key, expected, 'contraction')
    _bindings(core, base, packet.get('input_bindings'))
    good, remaining = core.validate_partition(packet.get('original_good_modes'),
                                              packet.get('original_remaining_modes'),
                                              1280, 1263, 18)
    matrix = _decode_upper(core, packet.get('inverse_free_border133_upper'), BORDER, BITS,
                           'inverse_free_border133_upper')
    _, reduced, paid, boundary = _reduce_fixed(core, base, matrix, FAR, BOUNDARY)
    return {
        'schema': REDUCTION_SCHEMA,
        'status': 'complete_fixed_reduction_conditional_on_inherited_source_contract',
        'eta': ['1', '1'], 'arithmetic_bits': BITS,
        'original_dimension': BORDER, 'reduced_dimension': FAR,
        'original_retained_indices': list(range(FAR)),
        'paid_original_indices': [131, 132], 'boundary_original_indices': list(range(112, 131)),
        'eliminated_inertia': [0, 21, 0],
        'paid_pivot_transcript': paid, 'boundary_pivot_transcript': boundary,
        'symmetric_encoding': ENCODING, 'reduced112_upper': _export_upper(reduced),
        'reduced_inertia': None, 'border_inertia': None,
        'original_good_modes': good, 'original_remaining_modes': remaining,
        'input_bindings': dict(packet['input_bindings']), 'padding_account': dict(PADDING),
        'reduction_bindings': {
            'contraction_result_sha256': expected_contraction_sha256,
            'reduction_source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        'premises': [
            'One original coherent realization, saved PX=B^T and G=BX>0 are inherited.',
            'Accepted W=C-epsilon I>0 and original Kc=-S^-1<0 are inherited source premises.',
            'Every removed paid/boundary pivot above was enclosed strictly negative.',
            'In(border133)=(0,21,0)+In(reduced112); retained inertia has not been counted.',
        ],
        'restricted_domination_claim': False, 'eta_budget_membership_claim': False,
        'tau_or_z_upper_claim': False, 'MAIN_or_RH_claim': False,
        'compiled_or_official_acceptance': False,
    }


def certify_reduced(core, base, reduction_raw, expected_reduction_sha256, plan_raw):
    """Complete one supplied covering112 plan, or raise with no final count.

    Every accepted pivot is invertible. An uncertain/zero-containing interval
    is failure of this certificate, never a singularity or nullity certificate.
    """
    _sources(core, base)
    packet = _authenticate(core, reduction_raw, expected_reduction_sha256, 'reduction result')
    for key, expected in {
        'schema': REDUCTION_SCHEMA,
        'status': 'complete_fixed_reduction_conditional_on_inherited_source_contract',
        'eta': ['1', '1'], 'arithmetic_bits': BITS,
        'original_dimension': BORDER, 'reduced_dimension': FAR,
        'original_retained_indices': list(range(FAR)),
        'paid_original_indices': [131, 132], 'boundary_original_indices': list(range(112, 131)),
        'eliminated_inertia': [0, 21, 0], 'symmetric_encoding': ENCODING,
        'padding_account': PADDING,
    }.items():
        core.equal_field(packet, key, expected, 'reduction')
    _bindings(core, base, packet.get('input_bindings'))
    reduction_bindings = packet.get('reduction_bindings')
    core.require(type(reduction_bindings) is dict, 'binding_mismatch', 'reduction_bindings')
    core.equal_field(reduction_bindings, 'reduction_source_sha256',
                     hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'reduction_bindings')
    core.require(type(plan_raw) is bytes, 'invalid_input_bytes', 'explicit immutable plan bytes')
    plan = core.parse_json(plan_raw)
    for key, expected in {
        'schema': PLAN_SCHEMA, 'purpose': 'inertia_research_only',
        'reduction_result_sha256': expected_reduction_sha256,
        'arithmetic_bits': BITS, 'original_retained_indices': list(range(FAR)),
    }.items():
        core.equal_field(plan, key, expected, 'plan')
    blocks = plan.get('blocks')
    _check_blocks(core, list(range(FAR)), blocks, complete=True)
    matrix = _decode_upper(core, packet.get('reduced112_upper'), FAR, BITS, 'reduced112_upper')
    _, active, inertia, transcript = _eliminate(core, base, matrix, list(range(FAR)), blocks)
    core.require(not active and sum(inertia) == FAR, 'incomplete_coverage', 'complete112 count')
    # This lower bound follows from the inherited bordered identity. A mismatch
    # is a contract inconsistency, not an admissible negative dimension.
    core.require(inertia[1] >= 56, 'inherited_contract_inconsistency',
                 'same-source restricted negative index must be nonnegative')
    restricted = [1205 + inertia[0], inertia[1] - 56, 0]
    return {
        'schema': RESULT_SCHEMA,
        'status': 'complete_signed_inertia_conditional_on_inherited_source_contract',
        'eta': ['1', '1'], 'arithmetic_bits': BITS,
        'original_retained_indices': list(range(FAR)), 'reduced_inertia': inertia,
        'border_inertia': [inertia[0], inertia[1] + 21, 0],
        'same_source_restricted_inertia': restricted,
        'same_source_restriction_positive_definite': inertia == [56, 56, 0],
        'pivot_transcript': transcript,
        'input_bindings': dict(packet['input_bindings']),
        'reduction_bindings': dict(reduction_bindings),
        'certificate_bindings': {
            'reduction_result_sha256': expected_reduction_sha256,
            'small_plan_sha256': hashlib.sha256(plan_raw).hexdigest(),
            'small_inertia_source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        'padding_account': dict(PADDING),
        'premises': list(packet['premises']) + [
            'No nullity is inferred from an uncertain pivot; all112 coordinates were covered.',
            'A finite same-source restriction count does not provide the independent eta budget.',
        ],
        'restricted_domination_claim': False, 'eta_budget_membership_claim': False,
        'tau_or_z_upper_claim': False, 'MAIN_or_RH_claim': False,
        'compiled_or_official_acceptance': False,
    }
