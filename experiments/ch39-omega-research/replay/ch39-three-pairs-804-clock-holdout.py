"""锁定单钟读法后，对三个不重叠的相邻零点对作一次批量检验。"""
from argparse import ArgumentParser
from pathlib import Path
from runpy import run_path


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    cli = ArgumentParser(description=__doc__)
    cli.add_argument("stage", choices=("forecast", "actual"))
    stage = cli.parse_args().stage
    evaluate = run_path(str(root / "ch39-tenth-804-clock-holdout.py"))["evaluate"]
    for source in (10, 12, 14):
        print("第几壳", source + 1, flush=True)
        evaluate(root, stage, previous_shell=source)
