"""三个等权子壳类组成的有理正混合，在给定的两组实际编号中共同选中中央。"""
from pathlib import Path
from runpy import run_path
import sys

# 键依次给出上下三站轨道数和极点标志、交叠数、上/下边界及概率分子。
SUPPORT = (
    ((1,True,1,False,0), "D","I", 14),
    ((1,False,2,True,0), "I","D", 79),
    ((3,True,6,False,3), "D","D", 1),
)


def verify(replay, tile=None):
    rank = run_path(str(replay / "ch39-exchangeable-finite-zero-rank.py"))
    rank["verify"](replay, tile, support=SUPPORT, denominator=94,
                   holdouts=((-4,36298,4813),(2,33560,12280)),
                   closest=(-4,31681), gap_bounds=("6017.03","6017.04"))


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent,
           int(sys.argv[2]) if len(sys.argv)>2 else None)
