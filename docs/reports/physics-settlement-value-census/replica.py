#!/usr/bin/env python3
"""Preview of the typed-source checks and counts in the main repository's
tools/lean-inspector/LeanInformationAuditRegAnalysis/AuricFib/Generic.lean, plus a
task-on-fiber check that the analyzer does not compute. Not the canonical analyzer.

Usage: python3 replica.py REQUEST_OR_BATCH.json
"""
import json
import sys


def validate(plan, finite):
    """The same input checks as Generic.analyze."""
    roles, n, kernels, layers = plan["roles"], finite["cardinality"], finite["kernels"], finite["layers"]
    initial, additions, task = plan["initial"], plan["additions"], plan["task"]

    def selected(k):
        return initial + [i for addition in additions[:k] for i in addition]

    assert task < len(roles) and all(i < len(roles) for i in selected(len(additions))), "fib.plan_role_scope"
    assert (len(finite["states"]) == n and len(finite["outputs"]) == len(roles)
            and all(len(o) == n for o in finite["outputs"]) and len(finite["task_outputs"]) == n
            and len(kernels) == len(roles) and len(layers) == len(additions) + 1), "fib.finite_presentation_shape"
    for table in kernels + layers:
        assert len(table) == n and all(len(row) == n for row in table), "fib.table_shape"
        assert all(table[i][j] == (table[i] == table[j]) for i in range(n) for j in range(n)), \
            "fib.table_kernel_reflection"
    for k, layer in enumerate(layers):
        assert all(layer[i][j] == all(kernels[r][i][j] for r in selected(k))
                   for i in range(n) for j in range(n)), "fib.layer_source_reflection"
    return selected


def report(request):
    plan, finite = request["plan"], request["finite"]
    selected = validate(plan, finite)
    roles, n, kernels, layers = plan["roles"], finite["cardinality"], finite["kernels"], finite["layers"]
    additions, task = plan["additions"], plan["task"]
    pairs = [(i, j) for i in range(n) for j in range(n) if i != j]
    total = len(pairs)

    def escape(readouts):
        return [(i, j) for i, j in pairs if all(kernels[r][i][j] for r in readouts)]

    def name(r):
        return roles[r]["name"]

    readouts = sorted(set(selected(len(additions))))
    full = escape(readouts)
    print(f"states={n}  ordered pairs={total}  readouts={[name(r) for r in readouts]}  task={name(task)}")
    print(f"escape pairs under all readouts: {len(full)}/{total}")
    print("unique capture / leave-one-out escape:")
    for r in readouts:
        without = escape([x for x in readouts if x != r])
        unique = [q for q in without if not kernels[r][q[0]][q[1]]]
        print(f"  {name(r):11s} unique={len(unique):4d}  leave-one-out={len(without)}")
    print("layered spectrum (" + " -> ".join(name(a[0]) for a in additions) + "):")
    remaining_before = total
    for k, layer in enumerate(layers):
        remaining = sum(layer[i][j] for i, j in pairs)
        captured = remaining_before - remaining
        collapsed = "  (collapsed)" if k and captured == 0 else ""
        print(f"  layer {k} {[name(r) for r in selected(k)]}: newly separated={captured}{collapsed}  remaining={remaining}")
        remaining_before = remaining
    mixed = [(i, j) for i, j in full if i < j and not kernels[task][i][j]]
    verdict = "sufficient" if not mixed else "not sufficient"
    print(f"task on unresolved fiber: unordered pairs with different task value = {len(mixed)} -> readouts {verdict}")
    fibers = {}
    for i in range(n):
        fibers.setdefault(tuple(finite["outputs"][r][i]["value"] for r in readouts), []).append(i)
    for key, members in sorted(fibers.items(), key=lambda item: -len(item[1])):
        values = sorted({finite["task_outputs"][i]["value"] for i in members})
        counts = [(v, sum(finite["task_outputs"][i]["value"] == v for i in members)) for v in values]
        flag = "  MIXED" if len(values) > 1 else ""
        prs = ",".join(str(finite["states"][i]["pr"]) for i in members)
        print(f"  fiber size {len(members):2d}{flag}: {' / '.join(key)}  PRs {prs}  task {counts}")


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: replica.py REQUEST_OR_BATCH.json")
    with open(sys.argv[1], encoding="utf-8") as handle:
        data = json.load(handle)
    for request in data if isinstance(data, list) else [data]:
        print("=== " + request["binding"].get("label", ""))
        report(request)


if __name__ == "__main__":
    main()
