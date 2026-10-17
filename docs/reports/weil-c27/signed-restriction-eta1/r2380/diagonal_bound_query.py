#!/usr/bin/env python3
"""One fixed-source scalar interval query; --self-test reads no real inputs.

All interval endpoints are integer multiples of 2^-640. Multiplication and
division round outward using integers. The last U upper is an exact Fraction.
The source-realization premises are inherited, not established by hashes.
No producer, solve, spectral routine, optional SPD clipping, or network call.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import io
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
import zlib


ROOT = Path(__file__).resolve().parent
RESEARCH = ROOT.parent
BITS = 640
SCALE = 1 << BITS
N, FAR, JET, BOUNDARY, RANK, WIDTH = 1280, 112, 56, 19, 131, 266
COMMIT = "f9610d3a5a732abaf68e733e50a7e34e66148ced"
TARGET_PATH = "docs/reports/weil-c27/joint-gamma-6765-rho32-v2-target.json.zlib"
TARGET_BLOB = "2245e47ffd0397a781885a82c6f48585bdd54802"
TARGET_SHA256 = "7546a915368b4a40e0e41ea0ea693403c239598c7a6e96b7b8b07e67e22e98c6"
TARGET_BYTES = 26662779
METRIC_BLOB = "10f33745b4d210073ea7392cdde212296399a5d4"
METRIC_SHA256 = "63bbcebbe00e69b64bf301654a7369ee8c7bb7ecf21dc589f9b408925f8f0d49"
METRIC_COMPRESSED_SHA256 = "9ba53b70a0f1cad25b9e92894eff93f1340dc31d1a9bf945c560b6dc6940d072"
METRIC_BYTES = 1012821
BASIS = "1/sqrt(d),sqrt(2/d)*cos(2*pi*n*x/d); no (-1)^n phase"
SOURCE_ID = "1ecb02691d9b99e4c181c3f8955ec84dd41328f277bb59bdffcf6d80412dd55e"
ORIGINAL_TARGET = "ea4a761f76b59a11b2d04faf80531901a5164f4145c2424a4e79889ae3f72e22"
BOUNDARY_ACTION = "94aec57dacacad6bae1698a7e643514e3dace976fddbcc0008b2b7fd5aab447c"
TARGET_PRODUCER = "10021a4faf837683e8430a7124241f3145e1bcfdce42aca3f817cdd4a491ea07"
METRIC_PRODUCER = "1166528a1418d2aae4d8026a6e5760b4de91b238b55966e138dda53d24437ae9"
INPUT_HASHES = {
    "arithmetic-scalars.json.zlib": "f993be1aa48824cbf22276ac2d80e834ca9e54e4519e2040f437ebb4f724990f",
    "middle-scalars.json.zlib": "55212d0181673c015de503928265509e594eb0f81b85e3afb34cc28a103d18db",
    "diagonal-floor-witness.json": "f45f3a333b50aa60c38c5b2ee9be416a165a1c719ba7dfb3ca6dfbe251bf1840",
    "far-coupling.json.zlib": "ac35edacd7f55de470586e3b8a5915bcaa49b9d2c4c9c2ee073eeaaa8c9c73e7",
    "inverse-residual-pilot.json.zlib": "6db5093bf78fe1dcf43cbbe1b39ef2b1faf4174cf7cf33fe796c0d10ad4d231b",
}
METRIC_PACKET_PINS = {
    Path(TARGET_PATH).name: TARGET_SHA256,
    "joint-gamma-6765-rho32-v2-pivots.json": "d441b9f9fc49338a121430a501043b76afcf7dffa9ec09c8d672aa9a2358b3f0",
    "full-gamma-unweighted-retained-row.json.zlib": "566cb68530f31dd6e26af076055bb368e2b2e82fdcc0857abb7b2b475a036f48",
    "full-gamma-finite-retained-row.json.zlib": "479b8da2bb1d490d0f4adf285a67797d733ec05d4dab33d88b6175946a89f28a",
}
DEFAULT_METRIC = RESEARCH / "published-metric-blobs-r2220/full-gamma-finite-good-block-metric.json"
DEFAULT_INTAKE = RESEARCH / "published-metric-payloads-verified-r2222.json"
INTEGER_TEXT = re.compile(r"-?(?:0|[1-9][0-9]*)\Z")


class BoundaryError(Exception):
    def __init__(self, code: str, detail: str):
        super().__init__(detail)
        self.code, self.detail = code, detail


def require(condition: bool, code: str, detail: str) -> None:
    if not condition:
        raise BoundaryError(code, detail)


def integer(value: object, where: str) -> int:
    require(type(value) is int, "invalid_integer", where)
    return value


def endpoint(value: object, where: str) -> int:
    require(isinstance(value, str) and INTEGER_TEXT.fullmatch(value) is not None,
            "invalid_endpoint", where)
    return int(value)


def ceil_div(a: int, b: int) -> int:
    require(type(b) is int and b > 0, "invalid_divisor", str(b))
    return -((-a) // b)


def ceil_fraction(x: Fraction) -> int:
    return ceil_div(x.numerator, x.denominator)


@dataclass(frozen=True)
class Interval:
    lo: int
    hi: int
    bits: int = BITS

    def __post_init__(self) -> None:
        require(type(self.lo) is int and type(self.hi) is int and self.lo <= self.hi,
                "invalid_interval", "ordered integer endpoints required")
        require(type(self.bits) is int and self.bits >= 0,
                "invalid_grid", "nonnegative integer grid required")

    @classmethod
    def decode(cls, pair: object, where: str, bits: int = BITS) -> Interval:
        require(type(pair) is list and len(pair) == 2, "invalid_pair", where)
        return cls(endpoint(pair[0], where + "[0]"), endpoint(pair[1], where + "[1]"), bits)

    @classmethod
    def rational(cls, lo: Fraction | int, hi: Fraction | int | None = None,
                 bits: int = BITS) -> Interval:
        lo = Fraction(lo)
        hi = lo if hi is None else Fraction(hi)
        require(lo <= hi, "invalid_interval", "reversed rational interval")
        scale = 1 << bits
        return cls((lo * scale).__floor__(), ceil_fraction(hi * scale), bits)

    def fractions(self) -> tuple[Fraction, Fraction]:
        scale = 1 << self.bits
        return Fraction(self.lo, scale), Fraction(self.hi, scale)

    def same_grid(self, other: Interval) -> None:
        require(self.bits == other.bits, "grid_mismatch", "interval operation")

    def __add__(self, other: Interval) -> Interval:
        self.same_grid(other)
        return Interval(self.lo + other.lo, self.hi + other.hi, self.bits)

    def __neg__(self) -> Interval:
        return Interval(-self.hi, -self.lo, self.bits)

    def __sub__(self, other: Interval) -> Interval:
        return self + (-other)

    def __mul__(self, other: Interval) -> Interval:
        self.same_grid(other)
        corners = (self.lo * other.lo, self.lo * other.hi,
                   self.hi * other.lo, self.hi * other.hi)
        scale = 1 << self.bits
        return Interval(min(corners) // scale, ceil_div(max(corners), scale), self.bits)

    def scale(self, k: int) -> Interval:
        integer(k, "integer interval scaling")
        return Interval(min(k * self.lo, k * self.hi), max(k * self.lo, k * self.hi), self.bits)

    def divide(self, n: int) -> Interval:
        require(type(n) is int and n > 0, "invalid_divisor", "positive exact integer required")
        return Interval(self.lo // n, ceil_div(self.hi, n), self.bits)

    def square(self) -> Interval:
        low = 0 if self.lo <= 0 <= self.hi else min(self.lo ** 2, self.hi ** 2)
        high = max(self.lo ** 2, self.hi ** 2)
        scale = 1 << self.bits
        return Interval(low // scale, ceil_div(high, scale), self.bits)

    def intersect(self, other: Interval) -> Interval:
        self.same_grid(other)
        lo, hi = max(self.lo, other.lo), min(self.hi, other.hi)
        require(lo <= hi, "empty_intersection", "source-real symmetric intervals disagree")
        return Interval(lo, hi, self.bits)

    def overlaps(self, other: Interval) -> bool:
        self.same_grid(other)
        return max(self.lo, other.lo) <= min(self.hi, other.hi)

    def is_zero(self) -> bool:
        return self.lo == self.hi == 0

    def closest_zero(self) -> Fraction:
        value = self.lo if self.lo > 0 else self.hi if self.hi < 0 else 0
        return Fraction(value, 1 << self.bits)

    def export(self) -> list[str]:
        return [str(self.lo), str(self.hi)]


def no_duplicates(pairs: list[tuple[str, object]]) -> dict:
    result = {}
    for key, value in pairs:
        require(key not in result, "duplicate_json_key", key)
        result[key] = value
    return result


def reject_constant(value: str) -> None:
    raise BoundaryError("invalid_json_constant", value)


def parse_json(raw: bytes) -> dict:
    obj = json.loads(raw.decode("utf-8"), object_pairs_hook=no_duplicates,
                     parse_float=Decimal, parse_constant=reject_constant)
    require(type(obj) is dict, "invalid_json_root", "object required")
    return obj


def equal_field(obj: dict, key: str, expected: object, where: str) -> None:
    require(key in obj and obj[key] == expected and type(obj[key]) is type(expected),
            "binding_mismatch", where + "." + key)


def decode_matrix(value: object, rows: int, cols: int, where: str) -> list[list[Interval]]:
    require(type(value) is list and len(value) == rows, "shape_mismatch", where)
    result = []
    for i, row in enumerate(value):
        require(type(row) is list and len(row) == cols, "shape_mismatch", f"{where}[{i}]")
        result.append([Interval.decode(pair, f"{where}[{i}][{j}]") for j, pair in enumerate(row)])
    return result


def validate_partition(positive: object, remaining: object, last: int,
                       positive_count: int, remaining_count: int) -> tuple[list[int], list[int]]:
    require(type(positive) is list and len(positive) == positive_count,
            "partition_mismatch", "positive_modes shape")
    require(type(remaining) is list and len(remaining) == remaining_count,
            "partition_mismatch", "remaining_modes shape")
    require(all(type(n) is int and 0 <= n <= last for n in positive + remaining),
            "partition_mismatch", "original labels must be integers in range")
    require(len(set(positive)) == len(positive) and len(set(remaining)) == len(remaining),
            "partition_mismatch", "duplicate original labels")
    require(not set(positive).intersection(remaining)
            and set(positive).union(remaining) == set(range(last + 1)),
            "partition_mismatch", "partition must be disjoint and complete")
    return list(positive), list(remaining)  # Deliberately preserve order.


def validate_metric(packet: dict) -> tuple[list[int], list[int]]:
    for key, expected in {
        "schema": "ACTUAL_C27_FINITE_GOOD_BLOCK_JOINT_METRIC_V1", "N": N, "M": 32768,
        "basis": BASIS, "source_bits": 384, "source_row": 0, "source_integer_sha256": SOURCE_ID,
        "grid_bits": BITS, "working_bits": 1024, "source_sha256": METRIC_PRODUCER,
        "joint_row_order": ["complete_unweighted_r", "finite_bJ"],
        "packet_sha256": METRIC_PACKET_PINS,
    }.items():
        equal_field(packet, key, expected, "metric")
    positive, remaining = validate_partition(packet.get("positive_modes"),
                                             packet.get("remaining_modes"), N, 1263, 18)
    decode_matrix(packet.get("inverse_actions"), 1263, 2, "metric.inverse_actions")
    decode_matrix(packet.get("joint_inverse_metric"), 2, 2, "metric.joint_inverse_metric")
    return positive, remaining


def validate_target(packet: dict) -> tuple[list[Interval], list[list[Interval]], list[list[Interval]], Interval]:
    for key, expected in {
        "schema": "ACTUAL_C27_STRUCTURED_TARGET_V1", "grid_bits": BITS,
        "producer_working_bits": 1024, "producer_source_sha256": TARGET_PRODUCER,
    }.items():
        equal_field(packet, key, expected, "target")
    meta = packet.get("metadata")
    require(type(meta) is dict, "shape_mismatch", "target.metadata")
    for key, expected in {
        "N": N, "M": 32768, "last_mode": 2048, "basis": BASIS,
        "far_feature_rank": FAR, "signed_update_rank": RANK, "total_generator_columns": WIDTH,
        "common_young_parameter": "1/32", "original_metadata_grid_bits": 512,
        "reused_far_allowance_grid_bits": BITS, "input_sha256": INPUT_HASHES,
        "original_target_sha256": ORIGINAL_TARGET, "joint_boundary_action_sha256": BOUNDARY_ACTION,
    }.items():
        equal_field(meta, key, expected, "target.metadata")
    bins = meta.get("positive_boundary_bins")
    require(type(bins) is list and len(bins) == BOUNDARY,
            "shape_mismatch", "nineteen boundary bins required")
    epsilon_numerator = (27852372618143874844 << 128) + (1 << 98)
    require(Interval.decode(meta.get("reused_far_allowance"), "reused_far_allowance")
            == Interval(epsilon_numerator, epsilon_numerator),
            "binding_mismatch", "complete transported epsilon")
    equal_field(meta, "complete_far_operator_allowance",
                ["27852372618143874843", "27852372618143874844"], "target.metadata")
    raw_diagonal = packet.get("diagonal")
    require(type(raw_diagonal) is list and len(raw_diagonal) == N + 1,
            "shape_mismatch", "target.diagonal")
    diagonal = [Interval.decode(pair, f"diagonal[{n}]") for n, pair in enumerate(raw_diagonal)]
    generators = decode_matrix(packet.get("generators"), N + 1, WIDTH, "generators")
    kernel = decode_matrix(packet.get("kernel"), RANK, RANK, "kernel")
    sigma = Interval.decode(packet.get("kappa"), "kappa")
    require(sigma.lo > 0, "binding_mismatch", "sigma must be positive")
    for n, row in enumerate(generators):
        for j in range(RANK):
            require(row[4 + j].scale(n * n).overlaps(row[4 + RANK + j]),
                    "generator_link_mismatch", f"original label {n}, feature {j}")
    return diagonal, generators, kernel, sigma


def build_plan(kernel: list[list[Interval]], far: int = FAR, jet: int = JET,
               boundary: int = BOUNDARY) -> tuple[list[tuple[int, int, Interval]], dict]:
    rank = far + boundary
    require(far == 2 * jet and jet > 0, "invalid_plan", "two far feature channels required")
    require(len(kernel) == rank and all(len(row) == rank for row in kernel),
            "shape_mismatch", "kernel plan")
    terms = []
    cross_zeros = mask_zeros = other_zeros = 0
    for i in range(rank):
        for j in range(i, rank):
            cross = i < far <= j
            masked = i < far and j < far and (i % jet) + (j % jet) >= jet
            if cross or masked:
                require(kernel[i][j].is_zero() and kernel[j][i].is_zero(),
                        "exact_zero_mismatch", f"kernel[{i},{j}] or transpose")
                cross_zeros += int(cross)
                mask_zeros += int(masked)
            shared = kernel[i][j].intersect(kernel[j][i])
            if shared.is_zero():
                other_zeros += int(not cross and not masked)
            else:
                terms.append((i, j, shared))
    return terms, {"cross_block_upper_zeros": cross_zeros,
                   "far_mask_upper_zeros": mask_zeros,
                   "other_exact_upper_zeros": other_zeros,
                   "quadratic_terms_per_row": len(terms)}


def quadratic(features: list[Interval], terms: list[tuple[int, int, Interval]]) -> Interval:
    require(bool(features), "shape_mismatch", "empty feature row")
    bits = features[0].bits
    require(all(x.bits == bits for x in features), "grid_mismatch", "feature row")
    lower = upper = 0
    for i, j, coefficient in terms:
        require(coefficient.bits == bits, "grid_mismatch", "quadratic coefficient")
        left, right = features[i], features[j]
        if i == j:
            product_lo = 0 if left.lo <= 0 <= left.hi else min(left.lo ** 2, left.hi ** 2)
            product_hi = max(left.lo ** 2, left.hi ** 2)
        else:
            products = (left.lo * right.lo, left.lo * right.hi,
                        left.hi * right.lo, left.hi * right.hi)
            product_lo, product_hi = min(products), max(products)
        corners = (coefficient.lo * product_lo, coefficient.lo * product_hi,
                   coefficient.hi * product_lo, coefficient.hi * product_hi)
        multiplicity = 1 if i == j else 2
        lower += multiplicity * min(corners)
        upper += multiplicity * max(corners)
    # All terms have common denominator 2^(3*bits); round just once.
    divisor = 1 << (2 * bits)
    return Interval(lower // divisor, ceil_div(upper, divisor), bits)


def reciprocal_pi(sigma: Interval) -> tuple[Fraction, Fraction]:
    lo, hi = sigma.fractions()
    lower, upper = max(lo / 2, Fraction(1, 4)), min(hi / 2, Fraction(1, 3))
    require(lower <= upper, "empty_pi_intersection", "sigma/2 intersect [1/4,1/3]")
    return lower, upper


def universal_upper(c_upper: Fraction, p: Interval, q: Interval, n: int,
                    t: tuple[Fraction, Fraction]) -> Fraction:
    require(type(n) is int and n > 0, "zero_or_invalid_label", "U requires original n>0")
    require(Fraction(1, 4) <= t[0] <= t[1] <= Fraction(1, 3),
            "invalid_pi_interval", "nearest-zero formula needs rational pi bounds")
    p0, q0 = p.closest_zero(), q.closest_zero()
    m = p0 + q0
    return c_upper - 2 * (abs(p0) + abs(q0)) + (t[1] if m >= 0 else t[0]) * m / n


def row_bound(position: int, n: int, ell: Interval, generators: list[Interval],
              terms: list[tuple[int, int, Interval]], t: tuple[Fraction, Fraction]) -> dict:
    c_diagonal = ell + quadratic(generators[4:4 + RANK], terms)
    result = {"original_good_list_position": position, "original_mode_label": n,
              "C_nn_interval": c_diagonal.export()}
    if n == 0:
        return dict(result, p_interval=None, q_interval=None, Uhat_upper_numerator=None,
                    zero_mode_marker="retained in full dimension; not divided by n or counted")
    a, b = generators[0], generators[1].divide(n)
    c = generators[0].divide(2) - generators[2]
    d = (generators[1].divide(2) - generators[3]).divide(n)
    p, q = a * d, c * b
    upper = universal_upper(c_diagonal.fractions()[1], p, q, n, t)
    return dict(result, p_interval=p.export(), q_interval=q.export(),
                Uhat_exact_rational=[str(upper.numerator), str(upper.denominator)],
                Uhat_upper_numerator=str(ceil_fraction(upper * (1 << ell.bits))))


def classify_rows(rows: list[dict], threshold: int = 78) -> tuple[int, str]:
    seen, count = set(), 0
    for position, row in enumerate(rows):
        n = integer(row.get("original_mode_label"), "row original label")
        require(n >= 0 and n not in seen, "row_label_mismatch", "duplicate or negative original label")
        require(row.get("original_good_list_position") == position,
                "row_order_mismatch", "saved original good-list order")
        seen.add(n)
        value = row.get("Uhat_upper_numerator")
        if n == 0:
            require(value is None, "zero_label_mismatch", "zero mode must be uncounted")
        else:
            count += int(endpoint(value, "Uhat upper numerator") <= 0)
    return count, "falsifier_certified" if count >= threshold else "not_certified"


def verify_sha256(raw: bytes, expected: str, where: str) -> str:
    digest = hashlib.sha256(raw).hexdigest()
    require(digest == expected, "sha256_mismatch", f"{where}: expected {expected}, observed {digest}")
    return digest


def git_blob_id(raw: bytes) -> str:
    digest = hashlib.sha1()
    digest.update(b"blob " + str(len(raw)).encode("ascii") + b"\0")
    digest.update(raw)
    return digest.hexdigest()


def git_read(repository: Path, args: list[str], audit: list[dict]) -> bytes:
    argv = ["git", "-c", "protocol.allow=never", "-C", str(repository), *args]
    env = dict(os.environ, GIT_NO_LAZY_FETCH="1", GIT_TERMINAL_PROMPT="0")
    try:
        proc = subprocess.run(argv, env=env, capture_output=True, timeout=60, check=False)
    except subprocess.TimeoutExpired as exc:
        audit.append({"argv": argv, "environment": {"GIT_NO_LAZY_FETCH": "1"},
                      "actual_exit": None, "timed_out": True,
                      "stdout_bytes": len(exc.stdout or b""),
                      "stderr": (exc.stderr or b"").decode("utf-8", errors="replace")})
        raise BoundaryError("git_read_timeout", str(args)) from exc
    except OSError as exc:
        audit.append({"argv": argv, "environment": {"GIT_NO_LAZY_FETCH": "1"},
                      "actual_exit": None, "launch_error": str(exc)})
        raise BoundaryError("git_launch_failure", str(args)) from exc
    audit.append({"argv": argv, "environment": {"GIT_NO_LAZY_FETCH": "1"},
                  "actual_exit": proc.returncode, "stdout_bytes": len(proc.stdout),
                  "stderr": proc.stderr.decode("utf-8", errors="replace")})
    require(proc.returncode == 0, "git_read_failed", f"{args}: exit {proc.returncode}")
    return proc.stdout


def read_fixed_target(repository: Path, audit: list[dict]) -> dict:
    commit = git_read(repository, ["rev-parse", "--verify", COMMIT + "^{commit}"], audit)
    require(commit == (COMMIT + "\n").encode(), "git_commit_mismatch", COMMIT)
    entry = git_read(repository, ["ls-tree", "-z", COMMIT, "--", TARGET_PATH], audit)
    expected = f"100644 blob {TARGET_BLOB}\t{TARGET_PATH}\0".encode()
    require(entry == expected, "git_path_mismatch", "fixed commit/path/blob identity")
    require(git_read(repository, ["cat-file", "-t", TARGET_BLOB], audit) == b"blob\n",
            "git_type_mismatch", TARGET_BLOB)
    require(git_read(repository, ["cat-file", "-s", TARGET_BLOB], audit)
            == f"{TARGET_BYTES}\n".encode(), "git_size_mismatch", TARGET_BLOB)
    raw = git_read(repository, ["cat-file", "blob", TARGET_BLOB], audit)
    require(len(raw) == TARGET_BYTES and git_blob_id(raw) == TARGET_BLOB,
            "git_blob_mismatch", "raw bytes and Git blob header SHA1")
    verify_sha256(raw, TARGET_SHA256, "compressed fixed target")
    decoder = zlib.decompressobj()
    decoded = decoder.decompress(raw) + decoder.flush()
    require(decoder.eof and not decoder.unused_data and not decoder.unconsumed_tail,
            "zlib_transport_mismatch", "one complete zlib stream required")
    return parse_json(decoded)


def read_authenticated_metric(path: Path, intake_path: Path, audit: list[dict]) -> dict:
    intake = parse_json(intake_path.read_bytes())
    equal_field(intake, "fixed_commit", COMMIT, "r2222 intake")
    packets = intake.get("packets")
    require(type(packets) is list, "intake_mismatch", "packet list")
    matches = [p for p in packets if type(p) is dict and p.get("name") == "full-gamma-finite-good-block-metric"]
    require(len(matches) == 1, "intake_mismatch", "one metric authentication entry")
    entry = matches[0]
    for key, expected in {
        "compressed_size": 303335, "sha256": METRIC_COMPRESSED_SHA256,
        "git_blob_sha1": METRIC_BLOB, "decoded_bytes": METRIC_BYTES,
        "decoded_sha256": METRIC_SHA256, "payload": str(DEFAULT_METRIC),
        "schema": str(RESEARCH / "published-metric-blobs-r2220/full-gamma-finite-good-block-metric-schema.json"),
    }.items():
        equal_field(entry, key, expected, "r2222 metric entry")
    raw = path.read_bytes()
    require(len(raw) == METRIC_BYTES, "metric_size_mismatch", str(path))
    digest = verify_sha256(raw, METRIC_SHA256, "existing decoded mode supplier")
    audit.append({"path": str(path), "bytes": len(raw), "verified_sha256": digest,
                  "authentication": "r2222 decoded supplier; no raw metric fetch or producer"})
    return parse_json(raw)


def run_query(args: argparse.Namespace, audit: list[dict]) -> dict:
    started = time.monotonic_ns()
    metric = read_authenticated_metric(args.metric, args.intake, audit)
    positive, remaining = validate_metric(metric)
    del metric
    packet = read_fixed_target(args.repository, audit)
    diagonal, generators, kernel, sigma = validate_target(packet)
    del packet
    terms, structure = build_plan(kernel)
    require(structure["cross_block_upper_zeros"] == 2128
            and structure["far_mask_upper_zeros"] == 3108 and len(terms) <= 3410,
            "kernel_structure_mismatch", "fixed 112+19 and far-mask counts")
    t = reciprocal_pi(sigma)
    rows = []
    for position, n in enumerate(positive):
        rows.append(row_bound(position, n, diagonal[n], generators[n], terms, t))
        if (position + 1) % 64 == 0:
            print(json.dumps({"stage": "scalar_rows", "completed": position + 1,
                              "total": len(positive), "elapsed_ns": time.monotonic_ns() - started}),
                  file=sys.stderr, flush=True)
    count, outcome = classify_rows(rows)
    return {
        "schema": "R2380_ORIGINAL_P_DIAGONAL_UNIVERSAL_FALSIFIER_V1",
        "status": outcome, "actual_process_exit": 0, "big_query_executed": True,
        "input_access": audit, "fixed_experiments_commit": COMMIT,
        "target_path": TARGET_PATH, "target_blob": TARGET_BLOB, "target_sha256": TARGET_SHA256,
        "metric_decoded_sha256": METRIC_SHA256, "source_integer_sha256": SOURCE_ID,
        "source_contract": {
            "basis": BASIS, "common_realization": "original shared symmetric kernel and generator rows",
            "sigma": "target kappa = 2/pi via original build/export and joint transport source chain",
            "sigma_source_blobs": ["09296f5808c8b078e812ce3445c71c26b4310767", "17b7ab717bc58da474575e00611a0d4caed9a74f"],
            "physical_far_scaling": "Ktarget_ff=(33/32)K; inherited binding, not newly proved by hashes",
            "negative_update_index_upper": "56 far directions + at most19 boundary directions =75",
            "negative_index_evidence": "inherited r2324 symmetric interval-family far-core certificate and original shared positive scaling",
            "paid_row_rank": 2, "boundary_SPD_needed_for_capacity_77": False,
            "source_binding_newly_proved": False, "node_feature_overlap_is_only_consistency": True,
        },
        "arithmetic": "integer dyadic intervals; quadratic accumulates on 2^-1920 then rounds once to640; final exact Fraction",
        "grid_bits": BITS, "sigma_interval": sigma.export(),
        "reciprocal_pi_interval": [[str(x.numerator), str(x.denominator)] for x in t],
        "epsilon": "27852372618143874844/2^512 + 2^-542; original full allowance retained",
        "saved_diagonal_role": "full-target ell; C_nn=ell+all131-feature quadratic",
        "boundary_SPD_clipping_enabled": False, "kernel_structure": structure,
        "positive_modes": positive, "remaining_modes": remaining, "zero_mode_is_good": 0 in positive,
        "rows": rows, "certified_nonpositive_upper_count": count, "threshold": 78,
        "conclusion": ("r2359 global-r,s diagonal scheme excluded for every eta>0; original P/stronger form unresolved"
                       if outcome == "falsifier_certified" else
                       "inconclusive only; no reverse implication or feasible eta/positivity established"),
        "eta_selected": False, "new_source_or_far_query": False, "inverse_or_solve": False,
        "producer_or_fixture_replay": False, "compiled_or_official_acceptance": False,
        "elapsed_ns": time.monotonic_ns() - started,
    }


def self_tests() -> dict:
    """Synthetic arithmetic only; never calls either real input reader."""
    completed = []

    def check(name: str, condition: bool) -> None:
        if not condition:
            raise AssertionError(name)
        completed.append(name)

    def rejects(name: str, call, code: str) -> None:
        try:
            call()
        except BoundaryError as exc:
            check(name, exc.code == code)
        else:
            raise AssertionError(name + ": did not reject")

    def iv(lo, hi=None, bits=6):
        return Interval.rational(Fraction(lo), None if hi is None else Fraction(hi), bits)

    quadrants = [((-3, -1), (2, 5), (-15, -2)),
                 ((-3, -1), (-4, -2), (2, 12)),
                 ((-2, 3), (-5, 7), (-15, 21)),
                 ((1, 3), (2, 4), (2, 12))]
    for index, (left, right, expected) in enumerate(quadrants):
        check(f"four_corner_product_{index}", (iv(*left) * iv(*right)) == iv(*expected))
    rational_product = iv(Fraction(1, 3), Fraction(2, 3)) * iv(Fraction(-5, 7), Fraction(3, 7))
    lo, hi = rational_product.fractions()
    check("rational_product_outward", lo <= Fraction(-10, 21) and hi >= Fraction(2, 7))
    check("square_crossing_zero", iv(-1, 2).square() == iv(0, 4))
    check("negative_ceil", ceil_div(-7, 3) == -2 and ceil_fraction(Fraction(-7, 3)) == -2)
    negative = Interval(-3, -1, 3).divide(3)
    check("negative_division_outward", negative == Interval(-1, 0, 3))
    check("symmetric_intersection", iv(1, 3).intersect(iv(2, 4)) == iv(2, 3))
    rejects("empty_symmetry_intersection", lambda: iv(1, 2).intersect(iv(3, 4)), "empty_intersection")
    rejects("grid_mismatch", lambda: iv(1).__add__(iv(1, bits=7)), "grid_mismatch")
    rejects("malformed_endpoint", lambda: Interval.decode(["1.0", "2"], "test"), "invalid_endpoint")
    rejects("reversed_endpoints", lambda: Interval.decode(["2", "1"], "test"), "invalid_interval")
    rejects("boolean_label", lambda: validate_partition([True], [0], 1, 1, 1), "partition_mismatch")
    check("original_order_preserved", validate_partition([3, 1, 4], [0, 2], 4, 3, 2)[0] == [3, 1, 4])
    rejects("duplicate_partition", lambda: validate_partition([1, 1, 4], [0, 2], 4, 3, 2), "partition_mismatch")
    rejects("overlapping_partition", lambda: validate_partition([1, 3, 4], [0, 1], 4, 3, 2), "partition_mismatch")
    rejects("missing_partition_label", lambda: validate_partition([1, 3], [0], 3, 2, 1), "partition_mismatch")
    rejects("matrix_shape", lambda: decode_matrix([], 1, 2, "synthetic"), "shape_mismatch")
    rejects("duplicate_json", lambda: parse_json(b'{"a":1,"a":2}'), "duplicate_json_key")
    rejects("nonfinite_json", lambda: parse_json(b'{"a":NaN}'), "invalid_json_constant")
    rejects("pin_mismatch", lambda: verify_sha256(b"synthetic", "0" * 64, "test"), "sha256_mismatch")
    check("git_blob_header_identity", git_blob_id(b"abc") == "f2ba8f84ab5c1bce84a7b441cb1959cfc7093b7f")

    zero = iv(0)
    small_kernel = [[zero for _ in range(5)] for _ in range(5)]
    small_kernel[0][0], small_kernel[4][4] = iv(2), iv(-4)
    small_kernel[0][2], small_kernel[2][0] = iv(1, 3), iv(2, 4)
    plan, stats = build_plan(small_kernel, far=4, jet=2, boundary=1)
    check("plan_uses_symmetric_intersection", (0, 2, iv(2, 3)) in plan)
    bad = [list(row) for row in small_kernel]
    bad[0][4] = iv(0, 1)
    rejects("crosszero_contains_zero_is_rejected", lambda: build_plan(bad, 4, 2, 1), "exact_zero_mismatch")
    bad = [list(row) for row in small_kernel]
    bad[1][1] = iv(-1, 1)
    rejects("mask_contains_zero_is_rejected", lambda: build_plan(bad, 4, 2, 1), "exact_zero_mismatch")
    bad = [list(row) for row in small_kernel]
    bad[0][2], bad[2][0] = iv(1), iv(2)
    rejects("plan_empty_symmetry", lambda: build_plan(bad, 4, 2, 1), "empty_intersection")

    features = [iv(1), iv(2), iv(3)]
    dense = [[2, 1, 0], [1, 3, 0], [0, 0, -4]]
    terms = [(i, j, iv(dense[i][j])) for i in range(3) for j in range(i, 3)]
    expected = sum(Fraction(dense[i][j]) * (i + 1) * (j + 1) for i in range(3) for j in range(3))
    check("contraction_matches_independent_dense", quadratic(features, terms) == iv(expected))
    uncertain_features = [iv(-1, 2), iv(-2, 1)]
    uncertain_coefficients = [iv(-2, 1), iv(-1, 3), iv(0, 2)]
    uncertain_terms = [(0, 0, uncertain_coefficients[0]), (0, 1, uncertain_coefficients[1]),
                       (1, 1, uncertain_coefficients[2])]
    enclosed = quadratic(uncertain_features, uncertain_terms).fractions()
    independent_values = []
    for x in (Fraction(-1), Fraction(0), Fraction(1, 2), Fraction(2)):
        for y in (Fraction(-2), Fraction(-1, 2), Fraction(0), Fraction(1)):
            for aa in (Fraction(-2), Fraction(1)):
                for bb in (Fraction(-1), Fraction(3)):
                    for cc in (Fraction(0), Fraction(2)):
                        independent_values.append(aa*x*x + bb*x*y + bb*y*x + cc*y*y)
    check("uncertain_contraction_independent_rationals", enclosed[0] <= min(independent_values)
          and max(independent_values) <= enclosed[1])
    tiny = Interval(1, 1, 3)
    check("quadratic_single_final_round", quadratic([tiny, tiny], [(0, 0, tiny), (1, 1, tiny)])
          == Interval(0, 1, 3))
    full_kernel = [[zero for _ in range(RANK)] for _ in range(RANK)]
    full_kernel[-1][-1] = iv(-7)
    full_plan, full_stats = build_plan(full_kernel)
    full_features = [zero] * RANK
    full_features[-1] = iv(1)
    check("all131_last_boundary_feature_retained", quadratic(full_features, full_plan) == iv(-7))
    check("fixed_exact_zero_counts", full_stats["cross_block_upper_zeros"] == 2128
          and full_stats["far_mask_upper_zeros"] == 3108)
    check("pi_source_interval", reciprocal_pi(iv(Fraction(1, 2), Fraction(2, 3)))
          == (Fraction(1, 4), Fraction(1, 3)))
    rejects("pi_empty_intersection", lambda: reciprocal_pi(iv(2)), "empty_pi_intersection")

    def samples(interval):
        lo, hi = interval.fractions()
        values = [lo + (hi - lo) * Fraction(k, 4) for k in range(5)]
        return values + ([Fraction(0)] if lo <= 0 <= hi else [])

    t = (Fraction(1, 4), Fraction(1, 3))
    cases = [(-2, -1, 3, 5), (1, 4, 2, 5), (-4, -2, -3, -1),
             (-3, 2, -5, -1), (-3, 4, -2, 5), (-2, -1, 1, 2)]
    for index, (pl, ph, ql, qh) in enumerate(cases):
        p, q = iv(pl, ph), iv(ql, qh)
        for n in (1, 2, 7):
            upper = universal_upper(Fraction(10), p, q, n, t)
            values = [Fraction(10) + tt * (pp + qq) / n - 2 * (abs(pp) + abs(qq))
                      for pp in samples(p) for qq in samples(q) for tt in (t[0], sum(t) / 2, t[1])]
            p0, q0 = p.closest_zero(), q.closest_zero()
            best_t = t[1] if p0 + q0 >= 0 else t[0]
            attained = Fraction(10) + best_t * (p0 + q0) / n - 2 * (abs(p0) + abs(q0))
            check(f"nearest_zero_box_max_{index}_n{n}", max(values) <= upper == attained)

    rows = [{"original_good_list_position": i, "original_mode_label": i + 1,
             "Uhat_upper_numerator": "0"} for i in range(78)]
    check("threshold_78_includes_equality", classify_rows(rows) == (78, "falsifier_certified"))
    check("threshold_77_is_inconclusive", classify_rows(rows[:-1]) == (77, "not_certified"))
    zero_row = {"original_good_list_position": 78, "original_mode_label": 0, "Uhat_upper_numerator": None}
    check("zero_label_not_counted", classify_rows(rows + [zero_row]) == (78, "falsifier_certified"))
    check("tiny_negative_ceil_zero", ceil_fraction(Fraction(-1, 2 * SCALE) * SCALE) == 0)
    check("tiny_positive_ceil_one", ceil_fraction(Fraction(1, 2 * SCALE) * SCALE) == 1)
    zrow = row_bound(0, 0, iv(3), [zero] * WIDTH, [], t)
    check("zero_row_no_division", zrow["Uhat_upper_numerator"] is None
          and zrow["C_nn_interval"] == iv(3).export())
    rejects("zero_label_direct_upper_rejected", lambda: universal_upper(Fraction(0), zero, zero, 0, t),
            "zero_or_invalid_label")
    error_stream = io.StringIO()
    exit_code = main(["--run-query"], stderr=error_stream)
    check("technical_failure_exit2_before_input_reads", exit_code == 2
          and json.loads(error_stream.getvalue())["status"] == "boundary_failure"
          and json.loads(error_stream.getvalue())["code"] == "cli_arguments"
          and json.loads(error_stream.getvalue())["input_access"] == [])
    return {"schema": "R2380_SYNTHETIC_SELF_TEST_V1", "status": "passed",
            "test_count": len(completed), "tests": completed,
            "real_input_readers_called": False, "big_query_executed": False,
            "source_or_mode_payload_read": False, "arithmetic": "integers/Fraction only"}


class Parser(argparse.ArgumentParser):
    def error(self, message: str) -> None:
        raise BoundaryError("cli_arguments", message)


def write_new(path: Path, result: dict) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


def main(argv: list[str] | None = None, stderr=None) -> int:
    stderr = sys.stderr if stderr is None else stderr
    audit = []
    args = None
    try:
        parser = Parser(description=__doc__)
        mode = parser.add_mutually_exclusive_group(required=True)
        mode.add_argument("--self-test", action="store_true")
        mode.add_argument("--run-query", action="store_true")
        parser.add_argument("--out", type=Path, help="fresh output JSON; required for --run-query")
        parser.add_argument("--repository", type=Path, default=Path("/Users/macstudio/omega/trureturing-experiments"))
        parser.add_argument("--metric", type=Path, default=DEFAULT_METRIC)
        parser.add_argument("--intake", type=Path, default=DEFAULT_INTAKE)
        args = parser.parse_args(argv)
        require(not args.run_query or args.out is not None, "cli_arguments", "--run-query requires --out")
        require(args.out is None or not args.out.exists(), "output_exists", "fresh output required")
        if args.self_test:
            saved_readers = read_fixed_target, read_authenticated_metric
            attempted_reads = []
            def blocked_reader(*unused_args, **unused_kwargs):
                attempted_reads.append(True)
                raise BoundaryError("self_test_input_guard", "real input reader prohibited in self-test")
            try:
                globals()["read_fixed_target"] = blocked_reader
                globals()["read_authenticated_metric"] = blocked_reader
                result = self_tests()
                require(not attempted_reads, "self_test_input_guard", "reader call attempted")
                result["real_input_reader_guard_enabled"] = True
            finally:
                globals()["read_fixed_target"], globals()["read_authenticated_metric"] = saved_readers
        else:
            result = run_query(args, audit)
        if args.out is not None:
            write_new(args.out, result)
        else:
            print(json.dumps(result, indent=2))
        return 0
    except Exception as exc:
        result = {"schema": "R2380_BOUNDARY_FAILURE_V1", "status": "boundary_failure",
                  "actual_process_exit": 2, "code": getattr(exc, "code", type(exc).__name__),
                  "detail": getattr(exc, "detail", str(exc)), "input_access": audit,
                  "mathematical_outcome_emitted": False}
        if args is not None and args.out is not None and not args.out.exists():
            try:
                write_new(args.out, result)
            except OSError as write_error:
                result["diagnostic_output_error"] = str(write_error)
        print(json.dumps(result, ensure_ascii=False), file=stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
