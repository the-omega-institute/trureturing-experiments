"""由单一跃迁概率或有符号 J 在最小割三窗内定位六面坐标。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


def verify(root):
    cp = run_path(str(root / "ch39-cp-ensemble-test.py"))
    ctx.prec, ctx.cap = 384,3
    data = cp["two"]["spectral_data"]()
    flavor = cp["flavor"]
    polynomials = cp["cp_data"](data,flavor["reference_matrix"]())
    reference = flavor["singlet"](11661)[1]
    amplitude = [flavor["singlet_at"]({724:z, 804:acb(1)})[1] for z in (acb(1),acb(0,1))]
    a,b = amplitude
    left,right = arb(".96619523"),arb(".99881043")
    assert a*arb(".256255").cos()-abs(b)*arb(".256255").sin() > left
    assert abs(acb(a,b)) < right
    assert left < reference < right
    u0 = cp["two"]["probabilities"](reference,data)[0][0]
    j0 = cp["cp_value"](reference,polynomials)
    assert u0 > 0 and j0 > 0
    pending = [(left,right)]
    leaves = 0
    while pending:
        low,high = pending.pop()
        box = low.union(high)
        up = -cp["two"]["probabilities"](box,data)[0][1]/u0
        cp_slope = -cp["cp_value"](box,polynomials,True)[1]/j0
        if arb("65.53") < up < arb("65.54") and arb("23.4") < cp_slope < 132:
            leaves += 1
            continue
        assert high-low>arb("1e-9")
        middle=(low+high)/2
        pending.extend(((low,middle),(middle,high)))
    print("六面坐标区间、导数分段数",left.str(12),right.str(12),leaves)
    print("两种归一化负导数的严格界", "65.53..65.54", "23.4..132")


if __name__=="__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
