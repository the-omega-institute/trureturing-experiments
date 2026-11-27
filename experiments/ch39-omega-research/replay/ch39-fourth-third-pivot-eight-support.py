"""八层低截止的全域动能条件只留下一个可能支撑。"""
from pathlib import Path
from runpy import run_path


def verify(root):
    core=run_path(str(root/"ch39-fourth-third-pivot-six-low-path-obstruction.py"))
    more,shifted=run_path(str(root/"ch39-fourth-third-pivot-nine-deformation.py"))["probe_list"]()
    allowed=[1,3,10,15,16,17,18,19]
    _,_,constants,high_error,bounds=core["integral_table"](
        root,more,shifted,(1334,2371,750000,17800000000,4200,5600,10000))
    core["integer_cover"](root,constants,high_error,bounds,8,except_rows=allowed)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
