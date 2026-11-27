"""右外壳上型对 92375 个九站支撑的轻质量排除。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    exceptions = [tuple(sorted((0,1,2,3,4,5,14) + pair)) for pair in combinations((11,12,18),2)]
    assert len(exceptions) == 3
    check = run_path(str(root/"ch39-up-eight-site-cover.py"))["verify"]
    check(root, 9, "ch39-up-nine-site-cover.tsv", ".00000629", exceptions)
