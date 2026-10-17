"""CLI for one prescribed fixed-source signed interval replay.

The caller supplies all local inputs and a complete pivot plan.
This source performs no Git/network access, plan selection or retry.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import sys
import time


CORE_SHA256 = "b278f8a358845ca436166b5595f8043173c7e09678d805d9d2c7e54d9d83d1c7"
REPLAY_SHA256 = "7806f5baab542868a8802cfe1119867657115f3c392183b3b893e37e3460ea98"


def load_pinned(path, expected, name):
    path = path.resolve()
    observed = hashlib.sha256(path.read_bytes()).hexdigest()
    if observed != expected:
        raise ValueError(f"source hash mismatch: {path}; observed {observed}")
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ValueError(f"cannot load pinned source: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core", type=Path, required=True)
    parser.add_argument("--target-compressed", type=Path, required=True)
    parser.add_argument("--metric-decoded", type=Path, required=True)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    # Avoid overwriting a previous receipt or any supplied input/source.
    protected = [args.core, args.target_compressed, args.metric_decoded, args.plan,
                 Path(__file__), Path(__file__).with_name("signed_interval_replay.py")]
    if args.output.resolve() in {p.resolve() for p in protected}:
        raise ValueError("output path aliases an input or source")
    if args.output.exists():
        raise FileExistsError(f"preserve existing receipt: {args.output}")
    started = time.monotonic_ns()
    core = None
    try:
        core = load_pinned(args.core, CORE_SHA256, "r2413_pinned_r2380")
        candidate = load_pinned(Path(__file__).with_name("signed_interval_replay.py"),
                                REPLAY_SHA256, "r2413_signed_replay")
        # Loading the plan precedes loading the target. Neither module selects pivots.
        plan_raw = args.plan.read_bytes()
        target_raw = args.target_compressed.read_bytes()
        metric_raw = args.metric_decoded.read_bytes()
        result = candidate.replay(core, target_raw, metric_raw, plan_raw)
        exit_code = 0
    except Exception as exc:
        known = core is not None and isinstance(exc, core.BoundaryError)
        result = {"schema": "R2413_REPLAY_FAILURE_V1", "status": "failed",
                  "failure_code": exc.code if known else type(exc).__name__,
                  "failure_detail": exc.detail if known else str(exc),
                  "inertia": None, "eta_budget_membership_claim": False,
                  "restricted_domination_claim": False, "tau_or_z_upper_claim": False,
                  "MAIN_or_RH_claim": False}
        exit_code = 2 if known else 1
    # ru_maxrss includes JSON decoding and the interpreter, in the native platform units.
    result["process_receipt"] = {
        "intended_process_exit": exit_code,
        "elapsed_ns": time.monotonic_ns() - started,
        "maximum_resident_set_size_native": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "maximum_resident_set_size_unit": "bytes" if sys.platform == "darwin" else "KiB on Linux; otherwise platform-defined",
        "platform": sys.platform,
        "execution_phase": "one prescribed inertia replay; no budget/MAIN conclusion"}
    # Exclusive create preserves earlier evidence; the parent records the actual shell exit.
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
