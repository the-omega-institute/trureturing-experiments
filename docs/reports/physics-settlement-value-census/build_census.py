#!/usr/bin/env python3
"""Build typed-source auricFibAnalysis requests for a census of research outcomes.

Inputs:
  --labels     labels.json: curated per-settlement readouts taken from PR bodies
  --main-repo  a git checkout of the-omega-institute/trureturing that contains every
               listed merge commit and its first parent (shallow clones: run
               `git fetch --depth=2 origin <sha>` for any commit reported missing)
  --out        output directory

Readouts computed here from the main repository:
  utility     strongest `utility:` header kind among the D5 modules added by the merge,
              read at the checkout's HEAD (certified-instance > checker > none)
  period      merge date bucket; exposure-time control
  downstream  task: some D5 module not added by the same merge imports an added module at HEAD

Outputs: physics-settlements-36.request.json (one request) and
physics-settlements-36.batch.json (two readout orders), both in the generic
typed-source format consumed by `tools/scripts/auric-fib-analysis.sh` in the main repository.
"""
import argparse
import collections
import json
import os
import re
import subprocess
import sys

PERIOD_SPLIT = "2026-10-09"
ROLES = [
    ("shape", "proof_shape of the settling theorem, as reported in the PR body (bind-only / content)"),
    ("admission", "module admission_basis reported in the PR body"),
    ("resolution", "Proved / Refuted / none (library result, no external question)"),
    ("utility", "strongest utility header kind among delivered D5 modules at HEAD"),
    ("track", "search track named in the PR body provenance (quick / deep / unlabelled)"),
    ("period", "merge date bucket; a control readout for exposure time"),
    ("downstream", "TASK: some D5 module outside this delivery imports a delivered module at HEAD (import is a proxy for use)"),
]
IMPORT = re.compile(r"^\s*(?:public\s+)?import\s+(D5\.[\w.]+)", re.M)
UTILITY = re.compile(r"utility:\s*(?:kind=([\w-]+)|(none))")


def fail(message):
    print(f"build_census: {message}", file=sys.stderr)
    sys.exit(2)


def git(repo, *args):
    result = subprocess.run(["git", "-C", repo, *args], capture_output=True, text=True)
    if result.returncode != 0:
        fail(f"git {' '.join(args)} failed: {result.stderr.strip()}")
    return result.stdout


def module_files(repo):
    files = {}
    for root, _, names in os.walk(os.path.join(repo, "D5")):
        for name in names:
            if name.endswith(".lean"):
                path = os.path.join(root, name)
                files[os.path.relpath(path, repo)[:-5].replace(os.sep, ".")] = path
    return files


def importers_of(files):
    importers = collections.defaultdict(set)
    for module, path in files.items():
        with open(path, encoding="utf-8", errors="ignore") as handle:
            for target in set(IMPORT.findall(handle.read(30000))):
                if target in files and target != module:
                    importers[target].add(module)
    return importers


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--labels", required=True)
    parser.add_argument("--main-repo", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    for path, kind in ((args.labels, "file"), (args.main_repo, "dir"), (args.out, "dir")):
        if (kind == "file" and not os.path.isfile(path)) or (kind == "dir" and not os.path.isdir(path)):
            fail(f"missing {kind}: {path}")
    labels = json.load(open(args.labels, encoding="utf-8"))
    settlements = labels.get("settlements")
    if not isinstance(settlements, list) or not settlements:
        fail("labels.json has no settlements")
    repo = os.path.abspath(args.main_repo)
    head = git(repo, "rev-parse", "HEAD").strip()
    for row in settlements:
        for sha in (row["merge_commit_sha"], row["merge_commit_sha"] + "^1"):
            if subprocess.run(["git", "-C", repo, "cat-file", "-e", sha + "^{commit}"],
                              capture_output=True).returncode != 0:
                fail(f"commit {sha} of PR {row['pr']} is absent; fetch it with --depth=2")
    files = module_files(repo)
    importers = importers_of(files)
    rows = []
    for row in settlements:
        sha = row["merge_commit_sha"]
        added = [line.split("\t")[1] for line in
                 git(repo, "diff", "--name-status", "--no-renames", sha + "^1", sha).splitlines()
                 if line.startswith("A\t") and line.split("\t")[1].startswith("D5/")
                 and line.endswith(".lean")]
        delivered = sorted({path[:-5].replace("/", ".") for path in added})
        present = [module for module in delivered if module in files]
        external = sorted(set().union(*(importers[m] for m in present)) - set(delivered)) if present else []
        kinds = set()
        for module in present:
            with open(files[module], encoding="utf-8") as handle:
                match = UTILITY.search(handle.read(4000))
            if match:
                kinds.add(match.group(1) or match.group(2))
        rows.append(dict(row,
                         utility=next((k for k in ("certified-instance", "checker", "none") if k in kinds), "unlabelled"),
                         period="early (2026-10-01..09)" if row["merged_at"] <= PERIOD_SPLIT else "late (2026-10-10..16)",
                         downstream="has external D5 importer" if external else "no external D5 importer",
                         delivered_modules=delivered, external_importers=external))
    fields = [name for name, _ in ROLES]
    n = len(rows)

    def table(keys):
        return [[all(rows[i][k] == rows[j][k] for k in keys) for j in range(n)] for i in range(n)]

    def request(order, label):
        return {
            "arena": "typed-source",
            "binding": {"kind": "research-outcome-census", "label": label,
                        "repository": "the-omega-institute/trureturing", "head": head,
                        "population": "merged physics settlement PRs created 2026-10-01..16 with proof_shape fields "
                                      "(36 of 39 merged; #12598 not physics; #17345, #17379 fields absent)",
                        "readout_sources": dict(ROLES),
                        "boundary": "labels are PR self-reports and repository reads; import is not live-path use; no source law"},
            "plan": {"available": True, "reason": "", "roles": [{"name": k, "meaning": d} for k, d in ROLES],
                     "initial": [], "additions": [[i] for i in order], "task": len(ROLES) - 1,
                     "scope": "36 research outcomes; equality kernels of the listed readouts"},
            "acquisition": {"kind": "finite",
                            "reason": "complete listed population of 36 merged PRs; labels extracted from PR bodies and the tree at head"},
            "finite": {"cardinality": n,
                       "states": [{"pr": r["pr"], "title": r["title"], "merged_at": r["merged_at"],
                                   "merge_sha": r["merge_commit_sha"], "delivered_modules": r["delivered_modules"],
                                   "external_importers": r["external_importers"]} for r in rows],
                       "outputs": [[{"value": r[k]} for r in rows] for k in fields],
                       "task_outputs": [{"value": r["downstream"]} for r in rows],
                       "kernels": [table([k]) for k in fields],
                       "layers": [table([fields[i] for i in order[:p]]) for p in range(len(order) + 1)]}}

    metrics_first = request([0, 1, 2, 3, 4, 5], "value readouts first, exposure control last")
    control_first = request([5, 0, 1, 2, 3, 4], "exposure control first")
    with open(os.path.join(args.out, "physics-settlements-36.request.json"), "w", encoding="utf-8") as handle:
        json.dump(metrics_first, handle, ensure_ascii=False, separators=(",", ":"))
        handle.write("\n")
    with open(os.path.join(args.out, "physics-settlements-36.batch.json"), "w", encoding="utf-8") as handle:
        json.dump([metrics_first, control_first], handle, ensure_ascii=False, separators=(",", ":"))
        handle.write("\n")
    print("| PR | merged | " + " | ".join(fields) + " | external importers |")
    print("|" + "---|" * (len(fields) + 3))
    for r in rows:
        print(f"| {r['pr']} | {r['merged_at']} | " + " | ".join(r[k] for k in fields) + f" | {len(r['external_importers'])} |")
    print(f"main repository HEAD {head}; {n} settlements", file=sys.stderr)


if __name__ == "__main__":
    main()
