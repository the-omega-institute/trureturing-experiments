"""给两只原始钟的部门间角差分别认证有限排除邻域。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def verify():
    engine = run_path(str(Path(__file__).with_name("ch39-shared-724-quark-gaps.py")))
    epsilon = arb(".005")
    for clock, endpoints, outer, middle in (
        (724, (".41","2.5"), ("2.65","13.9"), ("8.5","5.45")),
        (804, (".515",".95"), ("2.5","14.7"), ("11","4.15")),
    ):
        a, b = map(arb, endpoints)
        cases = [(arb(0),a,True),(a,b,False),(b,arb.pi(),True)]
        total = 0
        for sign in (1,-1):
            for left,right,is_outer in cases:
                u,d = outer if is_outer else middle
                for lo,hi in ((left-epsilon,left),(right,right+epsilon)):
                    for charge,index,limit in ((-4,int(is_outer),u),(2,int(not is_outer),d)):
                        total += engine["verify"](charge, sign, lo, hi, index, limit, shared_clock=clock)
        print("原始钟、扩展角宽、额外终格",clock,epsilon.str(6),total,flush=True)


if __name__ == "__main__":
    verify()
