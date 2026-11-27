"""在三个不重叠的一步留出上先更新两类边时钟，再揭晓完整第二枢轴指标。"""
from argparse import ArgumentParser
from pathlib import Path
from runpy import run_path


PAIRS = ((16, 17), (18, 19), (20, 21))


def batch(root, stage):
    score = run_path(str(root / "ch39-tenth-804-clock-holdout.py"))["score_bounds"]
    for start, target in PAIRS:
        assert target == start + 1
        print("起始、目标壳", start, target, flush=True)
        yield score(root, stage, previous_shell=start, update_724=True)


if __name__ == "__main__":
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("stage", choices=("forecast", "actual"))
    list(batch(Path(__file__).resolve().parent, cli.parse_args().stage))
