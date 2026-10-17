"""Candidate for one prescribed signed interval replay.

Revision r2429 extracts the identical state loop for independent small dense tests.

No CLI, default eta, pivot selection, scan, retry, solver, producer or Git IO.
The caller supplies the pinned r2380 module and already acquired byte inputs.
This returns In(eta C-epsilon I) only, conditional on the inherited source
contract. It never certifies budget membership, domination, tau, MAIN or RH.

The r2380 source is reused for interval arithmetic, source authentication,
the all131-feature diagonal recovery, and fixed target/metric validation.
Only signed reciprocals and the signed 1x1/2x2 Schur replay are added here.
Original r2413 preparation did not import or execute its candidate; r2429
validation results are recorded separately and do not certify real target input.
"""
from fractions import Fraction
import hashlib
from pathlib import Path
import zlib


CORE_SHA256 = "b278f8a358845ca436166b5595f8043173c7e09678d805d9d2c7e54d9d83d1c7"
SOURCE_PROOF_SHA256 = "b695d849fb1b1fa7047dab859051f8b54c167804f6322c20e3e4a0c67684c784"
PLAN_SCHEMA = "R2413_PRESCRIBED_SIGNED_REPLAY_PLAN_V1"
RESULT_SCHEMA = "R2413_SIGNED_FOUR_GENERATOR_INERTIA_V1"


def _rational(core, value, where):
    core.require(type(value) is list and len(value) == 2,
                 "invalid_rational", where)
    p = core.endpoint(value[0], where + ".numerator")
    q = core.endpoint(value[1], where + ".denominator")
    core.require(q > 0, "invalid_rational", where + ": positive denominator")
    return Fraction(p, q)


def _ceil_signed(a, b):
    # Valid for either sign of nonzero b, with Python's floor division.
    return -((-a) // b)


def _reciprocal(core, value, where):
    core.require(value.lo > 0 or value.hi < 0,
                 "uncertain_pivot", where + ": division interval contains zero")
    scale2 = 1 << (2 * value.bits)
    return core.Interval(scale2 // value.hi,
                         _ceil_signed(scale2, value.lo), value.bits)


def _scale_rational(core, value, scalar):
    numerators = (value.lo * scalar.numerator, value.hi * scalar.numerator)
    return core.Interval(min(numerators) // scalar.denominator,
                         core.ceil_div(max(numerators), scalar.denominator), value.bits)


def _divide_node_difference(core, value, difference):
    core.require(type(difference) is int and difference != 0,
                 "invalid_nodes", "distinct original squared labels required")
    return value.divide(difference) if difference > 0 else -value.divide(-difference)


def _lift(core, value, bits):
    core.require(bits >= value.bits, "invalid_grid", "no inward grid conversion")
    shift = bits - value.bits
    return core.Interval(value.lo << shift, value.hi << shift, bits)


def _skew_pair(core, left, right, eta_sigma):
    # left Omega_4 right^T, with Omega_4=eta*sigma*[[J,-J],[-J,0]].
    action = [eta_sigma * (right[1] - right[3]),
              -(eta_sigma * (right[0] - right[2])),
              -(eta_sigma * right[1]), eta_sigma * right[0]]
    total = core.Interval(0, 0, eta_sigma.bits)
    for a, b in zip(left, action):
        total = total + a * b
    return total


def _entry(core, i, j, generators, eta_sigma):
    numerator = _skew_pair(core, generators[i], generators[j], eta_sigma)
    return _divide_node_difference(core, numerator, i * i - j * j)


def _decode_bound_inputs(core, target_raw, metric_raw):
    core.require(type(target_raw) is bytes and type(metric_raw) is bytes,
                 "invalid_input_bytes", "caller must provide immutable bytes")
    core.require(len(target_raw) == core.TARGET_BYTES
                 and core.git_blob_id(target_raw) == core.TARGET_BLOB,
                 "git_blob_mismatch", "fixed compressed target")
    core.verify_sha256(target_raw, core.TARGET_SHA256, "fixed compressed target")
    core.require(len(metric_raw) == core.METRIC_BYTES,
                 "metric_size_mismatch", "existing decoded metric")
    core.verify_sha256(metric_raw, core.METRIC_SHA256, "existing decoded metric")
    decoder = zlib.decompressobj()
    decoded = decoder.decompress(target_raw) + decoder.flush()
    core.require(decoder.eof and not decoder.unused_data and not decoder.unconsumed_tail,
                 "zlib_transport_mismatch", "one complete target stream required")
    target = core.parse_json(decoded)
    metric = core.parse_json(metric_raw)
    labels, remaining = core.validate_metric(metric)
    ell, full_generators, kernel, sigma = core.validate_target(target)
    terms, structure = core.build_plan(kernel)
    core.require(structure["cross_block_upper_zeros"] == 2128
                 and structure["far_mask_upper_zeros"] == 3108 and len(terms) <= 3410,
                 "kernel_structure_mismatch", "same112+19 feature structure")
    return labels, remaining, ell, full_generators, sigma, terms


def _replay_state(core, labels, diagonal, generators, eta_sigma, blocks, dimension):
    """One shared signed Schur loop; dimensions are fixed by the authenticated caller."""
    bits = eta_sigma.bits
    active = list(labels)
    positive = negative = 0
    transcript = []
    for step, pivot_labels in enumerate(blocks):
        core.require(all(n in active for n in pivot_labels),
                     "invalid_plan", "a pivot label was already removed")
        rest = [n for n in active if n not in pivot_labels]
        pivot_rows = [list(generators[n]) for n in pivot_labels]
        if len(pivot_labels) == 1:
            k = pivot_labels[0]
            p = diagonal[k]
            inverse = [[_reciprocal(core, p, f"step {step}, label {k}")]]
            signs = (1, 0) if p.lo > 0 else (0, 1)
            pivot_matrix = [[p.export()]]
            determinant = p
        else:
            k, l = pivot_labels
            a, d = diagonal[k], diagonal[l]
            b = _entry(core, k, l, generators, eta_sigma)
            determinant = a * d - b.square()
            invdet = _reciprocal(core, determinant,
                                 f"step {step}, labels {k},{l}, determinant")
            inverse = [[d * invdet, -(b * invdet)],
                       [-(b * invdet), a * invdet]]
            if determinant.hi < 0:
                signs = (1, 1)
            else:
                trace = a + d
                core.require(trace.lo > 0 or trace.hi < 0,
                             "uncertain_pivot_sign",
                             f"step {step}, labels {k},{l}: positive determinant, uncertain trace")
                signs = (2, 0) if trace.lo > 0 else (0, 2)
            pivot_matrix = [[a.export(), b.export()], [b.export(), d.export()]]
        size = len(pivot_labels)
        for i in rest:
            # All c and pivot rows refer to the pre-update current Schur matrix.
            c = [_entry(core, i, k, generators, eta_sigma) for k in pivot_labels]
            multipliers = []
            for j in range(size):
                value = core.Interval(0, 0, bits)
                for t in range(size):
                    value = value + c[t] * inverse[t][j]
                multipliers.append(value)
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
        transcript.append({"step": step, "original_labels": list(pivot_labels),
                           "pivot_interval_matrix": pivot_matrix,
                           "determinant_interval": determinant.export(),
                           "inertia": [signs[0], signs[1], 0]})
        active = rest
    core.require(not active and positive + negative == dimension,
                 "incomplete_coverage", "no final count from a prefix")
    return positive, negative, transcript


def replay(core, target_raw, metric_raw, plan_raw):
    """Replay one supplied covering plan; return a complete inertia or failure.

    `core` must be the ordinary imported, pinned r2380 source module.
    This function does not acquire inputs, select eta or select pivot blocks.
    A successful result is research inertia only. An uncertain pivot aborts;
    it never reports singularity or a final negative count from a prefix.
    """
    core_path = Path(core.__file__)
    core.require(hashlib.sha256(core_path.read_bytes()).hexdigest() == CORE_SHA256,
                 "core_source_mismatch", "reuse the pinned r2380 source")
    core.require(type(plan_raw) is bytes, "invalid_input_bytes", "plan bytes required")
    plan = core.parse_json(plan_raw)
    for key, expected in {
        "schema": PLAN_SCHEMA, "purpose": "inertia_research_only",
        "target_sha256": core.TARGET_SHA256,
        "metric_decoded_sha256": core.METRIC_SHA256,
        "source_integer_sha256": core.SOURCE_ID,
        "source_structure_proof_sha256": SOURCE_PROOF_SHA256,
    }.items():
        core.equal_field(plan, key, expected, "plan")
    eta = _rational(core, plan.get("eta"), "plan.eta")
    core.require(eta > 0, "invalid_eta", "one explicit positive rational required")
    bits = core.integer(plan.get("arithmetic_bits"), "plan.arithmetic_bits")
    core.require(bits >= core.BITS, "invalid_grid", "at least the original640 grid")
    blocks = plan.get("blocks")
    core.require(type(blocks) is list and blocks, "invalid_plan", "pivot blocks required")
    flat = []
    for block in blocks:
        core.require(type(block) is list and len(block) in (1, 2)
                     and all(type(n) is int for n in block),
                     "invalid_plan", "blocks contain one or two original integer labels")
        flat.extend(block)
    core.require(len(flat) == 1263 and len(set(flat)) == 1263,
                 "incomplete_coverage", "exactly1263 distinct original labels required")

    labels, remaining, ell, full_generators, sigma, terms = _decode_bound_inputs(
        core, target_raw, metric_raw)
    core.require(set(flat) == set(labels), "partition_mismatch", "same saved good set")
    # The plan may permute pivots, but labels remain original; nodes are n^2.
    epsilon = core.Interval.rational(
        Fraction(27852372618143874844, 1 << 512) + Fraction(1, 1 << 542), bits=bits)
    eta_sigma = _scale_rational(core, _lift(core, sigma, bits), eta)
    diagonal, generators = {}, {}
    for n in labels:
        # Both old padding occurrences remain embedded in this exact recovery.
        c_diagonal = ell[n] + core.quadratic(full_generators[n][4:4 + core.RANK], terms)
        diagonal[n] = _scale_rational(core, _lift(core, c_diagonal, bits), eta) - epsilon
        generators[n] = [_lift(core, x, bits) for x in full_generators[n][:4]]
    # The large 266-column input rows are unnecessary during the replay itself.
    del ell, full_generators, terms
    positive, negative, transcript = _replay_state(
        core, labels, diagonal, generators, eta_sigma, blocks, 1263)
    return {
        "schema": RESULT_SCHEMA,
        "status": "complete_inertia_conditional_on_inherited_source_contract",
        "inertia": [positive, negative, 0],
        "eta": [str(eta.numerator), str(eta.denominator)],
        "arithmetic_bits": bits, "original_good_modes": labels,
        "original_remaining_modes": remaining, "pivot_transcript": transcript,
        "input_bindings": {"target_git_blob": core.TARGET_BLOB,
                           "target_sha256": core.TARGET_SHA256,
                           "metric_decoded_sha256": core.METRIC_SHA256,
                           "source_integer_sha256": core.SOURCE_ID,
                           "basis": core.BASIS, "core_source_sha256": CORE_SHA256,
                           "source_structure_proof_sha256": SOURCE_PROOF_SHA256,
                           "plan_sha256": hashlib.sha256(plan_raw).hexdigest()},
        "padding_account": {
            "C_initialization": "ell plus the complete shared131-feature quadratic",
            "embedded_charge": "a_orig + (1/32)*a_reuse; neither recovered separately nor altered",
            "additional_shift": "eta*C - epsilon*I; complete epsilon charged once",
            "epsilon_exact": "27852372618143874844/2^512 + 1/2^542"},
        "eta_budget_membership_claim": False,
        "restricted_domination_claim": False, "tau_or_z_upper_claim": False,
        "MAIN_or_RH_claim": False, "compiled_or_official_acceptance": False,
        "premises": [
            "Fixed hash authentication does not prove the inherited physical common-realization contract.",
            "C is the original coherent four-generator source with diagonal recovered from all131 features.",
            "Signed Schur displacement and congruence algebra are the ordinary r2401 derivations.",
            "No budget is needed to compute this standalone inertia; budget and final133-order actions/count remain separate."]}
