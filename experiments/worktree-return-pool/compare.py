#!/usr/bin/env python3
"""Compare canonical return/claim and remove/create on two authorized trees."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import threading
import time
import uuid


def command(argv, cwd):
    started = time.perf_counter()
    result = subprocess.run(list(map(str, argv)), cwd=cwd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True, timeout=3600)
    return dict(command=list(map(str, argv)), cwd=str(cwd), exit_code=result.returncode,
                wall_seconds=time.perf_counter() - started, stdout=result.stdout, stderr=result.stderr)


def checked(argv, cwd):
    result = command(argv, cwd)
    if result["exit_code"]:
        raise RuntimeError(json.dumps(result))
    return result["stdout"].strip()


def files(root):
    result = {}
    if root.exists():
        for parent, dirs, names in os.walk(root):
            dirs[:] = [name for name in dirs if name != ".git"]
            for name in names:
                path = Path(parent) / name
                stat = path.lstat()
                result[str(path.relative_to(root))] = (stat.st_size, stat.st_mtime_ns, stat.st_ino)
    return result


def disk(root):
    raw = checked(["df", "-k", root.parent], root.parent)
    return dict(df_k=raw, available_kib=int(raw.splitlines()[-1].split()[3]))


def changes(before, after):
    common = before.keys() & after.keys()
    modified = {name for name in common if before[name] != after[name]}
    return dict(files_created=len(after.keys() - before.keys()), files_removed=len(before.keys() - after.keys()),
                files_modified=len(modified),
                olean_files_written=sum(name.endswith(".olean") for name in (after.keys() - before.keys()) | modified))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--implementation", type=Path, required=True)
    parser.add_argument("--reuse-tree", type=Path, required=True)
    parser.add_argument("--recreate-tree", type=Path, required=True)
    parser.add_argument("--return-pr", type=int, required=True)
    parser.add_argument("--recreate-pr", type=int, required=True)
    parser.add_argument("--module", default="D5.S0.Carrier.TraceConjugation")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--skip-build-lock", action="store_true")
    parser.add_argument("--resume-builds", action="store_true")
    options = parser.parse_args()
    impl, reuse, recreate = [path.resolve() for path in (options.implementation, options.reuse_tree, options.recreate_tree)]
    if len({impl, reuse, recreate}) != 3 or any(str(path).startswith(("/tmp/", "/private/tmp/")) for path in (reuse, recreate)):
        raise ValueError("distinct non-temporary linked worktrees required")
    cli = ["dotnet", str(impl / "tools/StrataLint.Cli/bin/Release/net10.0/StrataLint.dll"), "worktree"]
    token = uuid.uuid4().hex[:12]
    report = (json.loads(options.output.read_text()) if options.resume_builds else
              dict(implementation_commit=checked(["git", "rev-parse", "HEAD"], impl),
                   module=options.module, phases=[], preparation=[], timestamp_unix=int(time.time())))
    options.output.parent.mkdir(parents=True, exist_ok=True)

    def save():
        staged = options.output.with_suffix(".tmp")
        staged.write_text(json.dumps(report, indent=2) + "\n")
        staged.replace(options.output)

    def phase(name, argv, target, cwd=None, setup=False):
        before, space = files(target), disk(target)
        result = command(argv, cwd or impl)
        after, final = files(target), disk(target)
        result.update(name=name, **changes(before, after), disk_before=space, disk_after=final,
                      available_kib_delta=final["available_kib"] - space["available_kib"],
                      built_lines=len(re.findall(r"\bBuilt\b", result["stdout"])),
                      cache_receipts=[line for line in result["stdout"].splitlines() if "LEAN_CACHE" in line])
        report["preparation" if setup else "phases"].append(result)
        save()
        print(json.dumps(dict(phase=name, exit_code=result["exit_code"], wall_seconds=result["wall_seconds"])), flush=True)
        if result["exit_code"]:
            raise RuntimeError(name + " failed; consult raw data")
        return result

    build_command = ["make", "lean", "LEAN_TARGETS=" + options.module]
    if options.skip_build_lock:
        build_command.append("LEAN_SKIP_LOCK=1")

    def conclude():
        report["final_heads"] = [checked(["git", "rev-parse", "HEAD"], path) for path in (reuse, recreate)]
        report["interpretation"] = dict(
            built_line_scope="Observed stdout count; canonical Lean door suppresses info-level Built messages.",
            file_scope="Entries under the subject, including ignored caches; Git directories and metadata outside the tree excluded. Modified means size/mtime/inode changed.",
            disk_scope="df available-space delta on the shared filesystem; concurrent host activity can contribute.",
            sample_scope="One pair and one concurrent exact-path claim trial, not a population timing estimate.")
        save()
        # Explicitly named cleanup is outside the comparison; every task ref remains.
        phase("cleanup-subjects", cli + ["remove", "--names", reuse.name + " " + recreate.name], reuse, setup=True)
        save()

    if options.resume_builds:
        if not options.skip_build_lock or report["module"] != options.module:
            raise ValueError("resumed build comparison requires the same module and explicit build-lock option")
        completed = {row["name"]: row for row in report["phases"]}
        for name in ("return", "reused-cache-door", "remove", "create", "recreated-cache-door"):
            if name not in completed or completed[name]["exit_code"]:
                raise ValueError("completed lifecycle and cache comparison required")
        clone = checked(["git", "rev-parse", "--path-format=absolute", "--git-common-dir"], impl)
        for path in (reuse, recreate):
            if checked(["git", "rev-parse", "--path-format=absolute", "--git-common-dir"], path) != clone:
                raise ValueError("foreign repository")
            if checked(["git", "rev-parse", "HEAD"], path) != report["source_commit"]:
                raise ValueError("subject revision changed")
            if (path / "lean-toolchain").read_text() != report["toolchain"]:
                raise ValueError("subject toolchain changed")
            if hashlib.sha256((path / "lake-manifest.json").read_bytes()).hexdigest() != report["manifest_sha256"]:
                raise ValueError("subject manifest changed")
        for row in report["phases"]:
            if row["name"] in ("reused-build", "recreated-build"):
                row["name"] = "default-lock-" + row["name"]
        report["build_options"] = ["LEAN_SKIP_LOCK=1"]
        save()
        phase("reused-build", build_command, reuse, cwd=reuse)
        phase("recreated-build", build_command, recreate, cwd=recreate)
        conclude()
        return

    report["build_options"] = ["LEAN_SKIP_LOCK=1"] if options.skip_build_lock else []
    common = checked(["git", "rev-parse", "--path-format=absolute", "--git-common-dir"], impl)
    for path, number in ((reuse, options.return_pr), (recreate, options.recreate_pr)):
        if checked(["git", "rev-parse", "--path-format=absolute", "--git-common-dir"], path) != common:
            raise ValueError("foreign repository")
        if checked(["git", "status", "--porcelain"], path):
            raise ValueError("ordinary work is present")
        pr = json.loads(checked(["gh", "pr", "view", str(number), "--json", "state,baseRefName,headRefName,headRefOid,mergeCommit"], path))
        branch = checked(["git", "branch", "--show-current"], path)
        if pr["state"] != "MERGED" or pr["baseRefName"] != "dev" or pr["headRefName"] != branch:
            raise ValueError("delivered subject identity required")
        report.setdefault("subjects", []).append(dict(path=str(path), branch=branch, pr=pr))
    listing = checked(["git", "worktree", "list", "--porcelain"], impl)
    if "branch refs/heads/lane/governance/idle-" in listing:
        raise ValueError("an unrelated idle candidate could enter the return sweep")
    source = checked(["git", "rev-parse", "HEAD"], reuse)
    report["source_commit"] = report["requested_base"] = source
    report["toolchain"] = (reuse / "lean-toolchain").read_text()
    report["manifest_sha256"] = hashlib.sha256((reuse / "lake-manifest.json").read_bytes()).hexdigest()
    save()
    # Only canonical removal/creation normalizes the second subject; its original
    # delivered task ref is retained. No reset, branch deletion or force is used.
    phase("normalize-remove", cli + ["remove", "--names", recreate.name], recreate, setup=True)
    phase("normalize-create", cli + ["--kind", "governance", "--name", "benchmark-initial-" + token,
          "--path", recreate, "--base", source], recreate, setup=True)
    for name, path in (("reuse", reuse), ("recreate", recreate)):
        phase("initial-cache-" + name, ["make", "lean-cache-ensure"], path, cwd=path, setup=True)
    report["same_condition"] = dict(heads=[checked(["git", "rev-parse", "HEAD"], path) for path in (reuse, recreate)],
        toolchains_equal=(reuse / "lean-toolchain").read_bytes() == (recreate / "lean-toolchain").read_bytes(),
        manifests_equal=(reuse / "lake-manifest.json").read_bytes() == (recreate / "lake-manifest.json").read_bytes(),
        cache_files=[sum(name.startswith(".lake/") for name in files(path)) for path in (reuse, recreate)])
    save()
    if report["same_condition"]["heads"] != [source, source]:
        raise ValueError("source commits differ")
    phase("return", cli + ["return", "--pr", options.return_pr, "--writers-joined"], reuse, cwd=reuse)
    before, space = files(reuse), disk(reuse)
    barrier = threading.Barrier(2)

    def contend(index):
        barrier.wait(timeout=60)
        return command(cli + ["--kind", "governance", "--name", "benchmark-claim-" + token + "-" + str(index),
                             "--path", reuse, "--base", source], impl)

    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(contend, range(2)))
    after, final = files(reuse), disk(reuse)
    report["concurrent_claim"] = dict(attempts=results, successes=sum(row["exit_code"] == 0 for row in results),
        **changes(before, after), disk_before=space, disk_after=final,
        available_kib_delta=final["available_kib"] - space["available_kib"])
    save()
    if report["concurrent_claim"]["successes"] != 1:
        raise RuntimeError("concurrent claim must have exactly one success")
    phase("reused-cache-door", ["make", "lean-cache-ensure"], reuse, cwd=reuse)
    phase("reused-build", build_command, reuse, cwd=reuse)
    phase("remove", cli + ["remove", "--names", recreate.name], recreate)
    phase("create", cli + ["--kind", "governance", "--name", "benchmark-recreated-" + token,
          "--path", recreate, "--base", source], recreate)
    phase("recreated-cache-door", ["make", "lean-cache-ensure"], recreate, cwd=recreate)
    phase("recreated-build", build_command, recreate, cwd=recreate)
    conclude()


if __name__ == "__main__":
    main()
