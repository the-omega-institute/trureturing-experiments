#!/usr/bin/env python3
"""Reconcile saved analyzer output, pair by pair, with a finite census.

This does not run or replace the Lean analyzer. It checks its saved output.
"""
import argparse
import contextlib
import hashlib
import io
import json
from pathlib import Path
import runpy


def digest(value):
    return hashlib.sha256(json.dumps(value, separators=(",", ":"), sort_keys=True).encode()).hexdigest()


def reconcile(request, actual):
    plan, finite = request["plan"], request["finite"]
    n, tables = finite["cardinality"], finite["kernels"]
    roles = [r["name"] for r in plan["roles"]]
    universe = {(i, j) for i in range(n) for j in range(n) if i != j}
    additions, initial, task = plan["additions"], plan["initial"], plan["task"]
    checks = []

    def check(name, expected, measured):
        checks.append(dict(name=name, expected=expected, actual=measured, equal=expected == measured))

    def pairs(name, expected, measured):
        given = [tuple(x) for x in measured]
        observed = set(given)
        checks.append(dict(name=name, expected_count=len(expected), actual_count=len(given),
                           expected_sha256=digest(sorted(expected)), actual_sha256=digest(sorted(given)),
                           missing=sorted(expected - observed), extra=sorted(observed - expected),
                           equal=expected == observed and len(given) == len(observed)))

    def selected(position):
        return initial + [r for a in additions[:position] for r in a]

    def escape(readouts):
        return {(i, j) for i, j in universe if all(tables[r][i][j] for r in readouts)}

    def rate(count):
        return dict(numerator=count, denominator=len(universe)) if universe else None

    all_roles = list(range(len(roles)))
    chosen = selected(len(additions))
    final, catalog_final = escape(chosen), escape(all_roles)
    conflicts = {(i, j) for i, j in final if not tables[task][i][j]}
    spectrum = actual["layered_spectrum"]
    sufficiency = actual["continuation_target"]["residual_sufficiency"]
    check("cardinality", n, actual["arena"]["cardinality"])
    check("ordered_pair_denominator", len(universe), spectrum["ordered_pair_denominator"])
    check("partition_count", len(universe), spectrum["partition_count"])
    check("source.supplied_binding", digest(request["binding"]), digest(actual["source_contract"]["supplied_binding"]))
    check("source.supplied_acquisition", request["acquisition"], actual["source_contract"]["supplied_acquisition"])
    check("source.binding_unverified", None, actual["source_contract"]["binding"])
    check("source.physical_law_certified", False, actual["source_contract"]["physical_law_certified"])
    check("plan.kernel_echo", digest(plan), digest(actual["kernel"]["structural"]))
    check("plan.spectrum_echo", digest(plan), digest(spectrum["structural"]))
    check("presentation_echo", finite, actual["arena"]["presentation"])
    # Avoid retaining another full input copy in this comparison summary.
    checks[-1] = {k: v for k, v in checks[-1].items() if k not in ("expected", "actual")}
    for r, table in enumerate(tables):
        check(f"kernel.{roles[r]}.outputs_reflect_table", True,
              all(table[i][j] == (finite["outputs"][r][i] == finite["outputs"][r][j])
                  for i in range(n) for j in range(n)))
    check("task_outputs_reflect_task_kernel", True,
          all(tables[task][i][j] == (finite["task_outputs"][i] == finite["task_outputs"][j])
              for i in range(n) for j in range(n)))
    check("kernel_tables_echo", digest(tables), digest(actual["kernel"]["finite_tables"]))
    check("kernel_outputs_echo", digest(finite["outputs"]), digest(actual["kernel"]["outputs"]))
    check("task_outputs_echo", digest(finite["task_outputs"]), digest(actual["continuation_target"]["outputs"]))
    check("catalog_full_relation", digest([[all(t[i][j] for t in tables) for j in range(n)] for i in range(n)]),
          digest(actual["kernel"]["full_relation"]))
    pairs("catalog_escape_pairs", catalog_final, actual["escape_pairs"])
    check("catalog_escape_rate", rate(len(catalog_final)), actual["escape_rate"])
    check("final_unresolved_count", len(final), spectrum["unresolved"]["count"])
    pairs("final_unresolved_pairs", final, spectrum["unresolved"]["pairs"])
    check("final_unresolved_rate", rate(len(final)), spectrum["unresolved"]["rate"])
    unique_rows = []
    check("unique_capture_role_count", len(all_roles), len(actual["unique_capture"]))
    for r, row in zip(all_roles, actual["unique_capture"]):
        without = escape([x for x in all_roles if x != r])
        unique = without - catalog_final
        check(f"unique.{roles[r]}.identity", r, row["readout"])
        check(f"unique.{roles[r]}.count", len(unique), row["count"])
        pairs(f"unique.{roles[r]}.pairs", unique, row["pairs"])
        check(f"unique.{roles[r]}.rate", rate(len(unique)), row["rate"])
        check(f"leave_one_out.{roles[r]}.count", len(without), row["leave_one_out_count"])
        pairs(f"leave_one_out.{roles[r]}.pairs", without, row["leave_one_out_pairs"])
        six_without = escape([x for x in chosen if x != r])
        unique_rows.append(dict(readout=roles[r], analyzer_catalog_unique=len(unique),
                                analyzer_catalog_leave_one_out=len(without),
                                preview_plan_unique=len(six_without - final) if r in chosen else None,
                                preview_plan_leave_one_out=len(six_without) if r in chosen else None))
    check("layer_count", len(additions) + 1, len(spectrum["layers"]))
    layer_rows = []
    for position, row in enumerate(spectrum["layers"]):
        remaining = escape(selected(position))
        before = universe if position == 0 else escape(selected(position - 1))
        captured = before - remaining
        check(f"layer.{position}.position", position, row["position"])
        check(f"layer.{position}.readouts", selected(position), row["readouts"])
        check(f"layer.{position}.count", len(captured), row["count"])
        pairs(f"layer.{position}.pairs", captured, row["pairs"])
        check(f"layer.{position}.rate", rate(len(captured)), row["rate"])
        check(f"layer.{position}.collapsed", position > 0 and remaining == before, row["collapsed"])
        layer_rows.append(dict(position=position, readouts=[roles[r] for r in selected(position)],
                               captured=len(captured), remaining=len(remaining), collapsed=row["collapsed"]))
    layer_escapes = [escape(selected(k)) for k in range(len(additions) + 1)]
    expected_first = {q: next((k for k, residual in enumerate(layer_escapes)
                              if q not in residual), None) for q in universe}
    observed_first = {tuple(x["pair"]): x["first_layer"] for x in spectrum["first_capture"]}
    check("first_capture.address_count", len(universe), len(spectrum["first_capture"]))
    check("first_capture.all_pair_addresses_and_nulls", True, expected_first == observed_first)
    determined = next((k for k in range(len(additions) + 1)
                       if all(tables[task][i][j] for i, j in escape(selected(k)))), None)
    check("residual.task", task, sufficiency["task"])
    check("residual.sufficient", not conflicts, sufficiency["sufficient"])
    check("residual.conflict_count", len(conflicts), sufficiency["conflict_count"])
    pairs("residual.conflict_pairs", conflicts, sufficiency["conflict_pairs"])
    check("residual.conflict_rate", rate(len(conflicts)), sufficiency["conflict_rate"])
    check("residual.determined_at_layer", determined, sufficiency["determined_at_layer"])
    prs = [x["pr"] for x in finite["states"]]
    return dict(label=request["binding"]["label"], passed=all(x["equal"] for x in checks),
                states=n, ordered_pairs=len(universe), final_unresolved=len(final),
                conflict_count=len(conflicts), sufficient=not conflicts, determined_at_layer=determined,
                conflict_pr_pairs=sorted([prs[i], prs[j]] for i, j in conflicts),
                scope_comparison=unique_rows, layers=layer_rows, checks=checks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batch", type=Path)
    parser.add_argument("analyzer_output", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--registered", action="store_true")
    args = parser.parse_args()
    requests, actual = json.loads(args.batch.read_text()), json.loads(args.analyzer_output.read_text())
    if len(requests) != len(actual):
        raise SystemExit("batch result cardinality differs")
    preview = runpy.run_path(str(Path(__file__).with_name("replica.py")))
    stream = io.StringIO()
    for request in requests:
        preview["validate"](request["plan"], request["finite"])
    if args.registered:
        with contextlib.redirect_stdout(stream):
            for request in requests:
                print("=== " + request["binding"].get("label", ""))
                preview["report"](request)
    report = dict(batch_sha256=hashlib.sha256(args.batch.read_bytes()).hexdigest(),
                  analyzer_output_sha256=hashlib.sha256(args.analyzer_output.read_bytes()).hexdigest(),
                  comparisons=[reconcile(r, a) for r, a in zip(requests, actual)],
                  scope_note="Top-level escape/unique/leave-one-out use all seven catalog roles, including task; "
                             "the preview and planned layers use only the six selected observation roles. "
                             "These are different definitions. Compare each against its own role set; "
                             "do not substitute 92 for the 102 planned final unresolved pairs.")
    if args.registered:
        report["stored_preview_exact_match"] = stream.getvalue() == Path(__file__).with_name("replica-preview.txt").read_text()
        expected_conflicts = sorted([[a, b] for a, b in [(15354, 16702), (16702, 16902),
                                 (16702, 17145), (16702, 17605), (16971, 17778)]
                                 for a, b in [(a, b), (b, a)]])
        report["registration_matches"] = report["batch_sha256"] == "8a26793329161a9b0654d5976e187c8e2d167b72cebd273fab57e51d6d90004d" and all(
            [x["states"], x["ordered_pairs"], x["final_unresolved"], x["conflict_count"],
             x["sufficient"], x["determined_at_layer"], x["conflict_pr_pairs"]] ==
            [36, 1260, 102, 10, False, None, expected_conflicts] for x in report["comparisons"])
    report["passed"] = all(x["passed"] for x in report["comparisons"]) and all(
        report.get(k, True) for k in ["stored_preview_exact_match", "registration_matches"])
    args.out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "comparisons"}, ensure_ascii=False))
    if not report["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
