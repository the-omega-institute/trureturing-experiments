#!/usr/bin/env python3
"""Replay the retained E2a readings and check their evidence boundaries.

This verifies published data and scoring, never solves the frozen event.
The pinned main checkout is read only; Python runs with -B.
"""

import argparse
from datetime import datetime
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--main-repo", type=Path, required=True)
    args = parser.parse_args()
    data_dir = Path(__file__).resolve().parent
    report = data_dir.parent
    root = report.parents[2]
    metadata = json.loads((data_dir / "data.json").read_text())
    sources = json.loads((data_dir / "sources.json").read_text())
    base = sources["frozen_input_commit"]
    main_repo = args.main_repo.resolve()
    env = dict(os.environ, GIT_OPTIONAL_LOCKS="0", PYTHONDONTWRITEBYTECODE="1")

    def git(cwd, *arguments):
        return subprocess.check_output(["git", *arguments], cwd=cwd, env=env)

    assert git(main_repo, "rev-parse", "HEAD").decode().strip() == metadata["main_baseline"]
    code_path = "tools/scripts/agent/research_readouts/readouts.py"
    blob = git(main_repo, "rev-parse", metadata["main_baseline"] + ":" + code_path).decode().strip()
    assert blob == "5e4d264cb223a2c8287e5472179fccfaca887c09"
    assert (main_repo / code_path).read_bytes() == git(main_repo, "show", metadata["main_baseline"] + ":" + code_path)

    hashes = {}
    for name, expected in metadata["sha256"].items():
        hashes[name] = hashlib.sha256((data_dir / name).read_bytes()).hexdigest()
        assert hashes[name] == expected, name
    online = (data_dir / "readouts-online.json").read_bytes()
    assert online == (data_dir / "readouts-offline.json").read_bytes()
    replay_env = dict(env, PATH=str(Path(sys.executable).resolve().parent) + os.pathsep + os.environ["PATH"])
    replay = subprocess.run(
        ["make", "-s", "-C", "tools", "research-readouts", "AS_OF=" + metadata["as_of"]],
        cwd=main_repo, env=replay_env, input=(data_dir / "issue-18101.rest.json").read_bytes(),
        capture_output=True, check=True,
    )
    assert not replay.stderr, replay.stderr
    assert replay.stdout == online
    assert metadata["online"]["exit_code"] == metadata["offline"]["exit_code"] == 0

    snapshot = json.loads((data_dir / "issue-18101.rest.json").read_text(), parse_float=Decimal)
    assert len(snapshot["issues"]) == 1 and snapshot["evidence"] == []
    item = snapshot["issues"][0]
    issue, comments = item["issue"], item["comments"]
    assert issue["number"] == 18101
    assert len(comments) == issue["comments"] == metadata["rest"]["retained_comment_count"] == 4
    assert sum(metadata["rest"]["page_sizes"]) == len(comments)
    by_id = {comment["id"]: comment for comment in comments}
    assert len(by_id) == len(comments)
    assert list(by_id) == metadata["rest"]["comment_ids"]

    def block(comment_id, kind):
        comment = by_id[comment_id]
        assert comment["created_at"] == comment["updated_at"]
        matches = re.findall(r"^```" + kind + r"\n(.*?)\n```", comment["body"], re.M | re.S)
        assert len(matches) == 1
        return json.loads(matches[0], parse_float=Decimal)

    prereg = block(6231250391, "research-prereg-v1")
    terminal = block(6231318725, "research-terminal-v1")
    decisions = json.loads((report / "decisions.json").read_text())
    frozen = decisions["slots"]["residual"]["event"]
    reading = json.loads(online)
    unit = reading["units"][0]
    assert prereg["event"] == terminal["event"] == unit["event"] == frozen
    assert all(seat["event"] == frozen for seat in prereg["predictions"])
    assert decisions["slots"]["boundary"]["status"] == "infeasible"
    assert decisions["slots"]["boundary"]["event"] is None
    assert prereg["scoring_rule"] == {"name": "binary-brier", "version": 1}
    assert Fraction(prereg["surprise_threshold"]) == Fraction(9, 16)

    def utc(value):
        return datetime.strptime(value, "%Y-%m-%dT%H:%M:%SZ")

    registered = utc(by_id[6231250391]["created_at"])
    started, settled = map(utc, [terminal["execution_started_at"], terminal["settled_at"]])
    deadline = utc(prereg["deadline"])
    assert utc(decisions["frozen_at"]) < registered < started < settled
    assert settled < utc(by_id[6231287063]["created_at"]) < utc(by_id[6231318594]["created_at"])
    assert utc(by_id[6231318594]["created_at"]) < utc(by_id[6231318725]["created_at"]) < deadline
    publication = re.search(r"preselected publication second ([0-9T:Z-]+)", by_id[6231250391]["body"])[1]
    assert (deadline - utc(publication)).total_seconds() == 3600
    assert (deadline - registered).total_seconds() == 3599
    assert (started - registered).total_seconds() == 36
    assert (settled - started).total_seconds() == 91
    assert terminal["status"] == "refuted" and terminal["shape_evidence"] is None
    evidence = terminal["evidence"]
    assert evidence["url"] == by_id[6231287063]["html_url"]
    assert evidence["excerpt"] in by_id[6231287063]["body"]

    # Recompute prediction scoring only, using the independently evaluated a=0.
    losses = {seat["seat"]: Fraction(seat["probability"]) ** 2 for seat in prereg["predictions"]}
    surprise = sum(losses.values(), Fraction()) / len(losses)
    mean = sum((Fraction(seat["probability"]) for seat in prereg["predictions"]), Fraction()) / len(losses)
    assert surprise == Fraction(401, 2000000) and mean == Fraction(21, 2000)
    assert all(Fraction(score["brier"]) == losses[score["seat"]] for score in unit["scores"])
    assert Fraction(unit["surprise"]) == surprise and Fraction(unit["aggregate_probability"]) == mean
    assert unit["y"] == int(surprise >= Fraction(9, 16)) == 0 and unit["z"] == 0
    assert unit["x"] is None and unit["mode"] is None
    assert unit["evidence_status"] == "missing-shape" and unit["score_status"] == "scored"
    assert unit["model_families"] == ["GPT-6"]
    evaluator = by_id[6231318594]["body"]
    assert all(value in evaluator for value in ["a=0", "1/1000000", "1/2500", "401/2000000", "21/2000", "非已编译声明的活路径"])
    batch = reading["batch"]
    assert batch["unit_count"] == batch["scored_units"] == 1
    assert batch["score_coverage"] == "1" and batch["mode_coverage"] == "0"
    assert batch["mode_denominator"] == 0 and all(value == 0 for value in batch["five_counts"].values())
    assert all(batch[key] is None for key in ["X", "Y", "Z"])
    assert batch["missing_counts"] == {"answer_units": 0, "mode_units": 1, "prediction_units": 0, "shape_units": 1}

    roles = sources["role_evidence"]
    flights = [seat["flight"] for seat in roles["predictors"]] + [roles["solver"]["flight"], roles["evaluator"]["flight"]]
    assert len(set(flights)) == 4
    assert roles["solver"]["saw_predictions"] is False
    assert "skills/autonomous-research/SKILL.md" in roles["solver"]["skills_loaded"]
    evaluator_fixed = utc(roles["evaluator"]["checks_fixed_before_reading_predictions_at"])
    evaluator_completed = utc(roles["evaluator"]["completed_at"])
    assert utc(by_id[6231287063]["created_at"]) < evaluator_fixed < evaluator_completed < utc(by_id[6231318594]["created_at"])
    evaluator_enclosure = int((utc(by_id[6231318594]["created_at"]) - utc(by_id[6231287063]["created_at"])).total_seconds())
    assert evaluator_enclosure == 227 and evaluator_enclosure < 600

    relative_report = report.relative_to(root).as_posix()
    existing = git(root, "ls-tree", "-r", "--name-only", base, "--", relative_report).decode().splitlines()
    protected = [path for path in existing if path != relative_report + "/README.md"]
    diff = subprocess.run(["git", "diff", "--exit-code", base, "--", *protected], cwd=root, env=env, capture_output=True, check=True)
    assert diff.stdout == b""
    old_readme = git(root, "show", base + ":" + relative_report + "/README.md")
    assert (report / "README.md").read_bytes().startswith(old_readme)
    changes = git(root, "diff", "--name-status", base).decode().splitlines()
    assert all(line.startswith("A\t" + relative_report + "/") or line == "M\t" + relative_report + "/README.md" for line in changes)
    print(json.dumps({
        "format": "e2a-result-verification-v1", "as_of": metadata["as_of"],
        "l2_blob": blob, "sha256": hashes, "offline_replay_exit_code": replay.returncode,
        "online_offline_replay_byte_equal": True, "complete_rest_comment_count": len(comments),
        "event_and_prediction_text_unchanged": True, "registered_window_seconds": 3599,
        "prediction_to_execution_seconds": int((started - registered).total_seconds()),
        "execution_to_settlement_seconds": int((settled - started).total_seconds()),
        "evaluator_enclosure_seconds": evaluator_enclosure,
        "terminal_before_deadline_seconds": int((deadline - utc(by_id[6231318725]["created_at"])).total_seconds()),
        "evaluator_arithmetic_matches_l2": True, "x_and_mode_unmeasured": True,
        "distinct_prediction_solver_evaluator_flights": 4,
        "protected_base": base, "protected_diff_exit_code": diff.returncode,
        "protected_paths": protected, "readme_is_append_only": True,
        "new_event_mathematics_executed": False,
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
