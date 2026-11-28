"""在两个钟的整个相位域，由单向乘积估计统一改进两种等能量条件的振幅下界。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb,ctx


def verify(root):
    run_path(str(root/"ch39-proton-equal-energy-torus.py"))["verify"]()
    run_path(str(root/"ch39-proton-induced-energy-exchange.py"))["verify"](root)
    ctx.prec=192
    for numerator,denominator,full_limit,three_limit in (
            (3,8,".00097817",".00148132"),(11,19,".00022594",".00067620")):
        dominant_cap=(arb(numerator)/denominator).sqrt()
        total_cap=1+2*dominant_cap-2*dominant_cap.sqrt()
        global_bound=arb(".0052")-arb(".0064")*total_cap
        three_bound=arb(".006")-arb(".00685")*total_cap
        assert global_bound>arb(full_limit) and three_bound>arb(three_limit)
        print("单向乘积上界、全相位及三枚零点的实振幅下界",dominant_cap.str(13),
              global_bound.str(13),three_bound.str(13))


if __name__=="__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
