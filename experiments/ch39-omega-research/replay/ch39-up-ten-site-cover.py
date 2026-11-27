"""十站轻上夸克检验中由逐层界排除的支撑。"""
from pathlib import Path
from runpy import run_path


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    lines = (root/"ch39-up-ten-site-exceptions.tsv").read_text().splitlines()
    exceptions = [tuple(map(int,row.split())) for row in lines]
    assert len(set(exceptions)) == 32
    assert all(len(s)==10 and s==tuple(sorted(set(s))) and min(s)>=0 and max(s)<19 for s in exceptions)
    check = run_path(str(root/"ch39-up-eight-site-cover.py"))["verify"]
    check(root, 10, "ch39-up-ten-site-cover.tsv", ".00000629", exceptions)
