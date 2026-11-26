"""对称且逐字概率随字长不增时，直接矩末圈份额的精确极值证书。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat, ctx


def verify(root):
    geometry = run_path(str(root/"ch39-mincut-real-triangle-tomography.py"))
    ctx.prec = 256
    u = geometry["unit"]((108469,109))
    theta = u.arg()
    rows = []
    for q in range(1,106):
        running = arb(1)
        row = [running]
        for r in range(1,107):
            running += 2*(q*r*theta).cos()
            row.append(running/(2*r+1))
        rows.append(row)
    cases = (
        (-1, (1,105,106), ".0637103937438668", ".0637103937438670"),
        (1, (0,87,89,91,93,95,97,99,101,103,105), ".0136255369474", ".0136255369477"),
    )
    for sign, active, lo, hi in cases:
        X = arb_mat([[row[r] for r in active] for row in rows])
        gram = X.transpose()*X
        assert not gram.det().contains(0)
        coefficients = gram.solve(arb_mat([[sign*rows[-1][r]] for r in active]))
        c = coefficients.entries()
        assert all(z > 0 for z in c)
        fitted = (X*coefficients).entries()
        residual = [arb(sign*int(q == 105))-v for q,v in enumerate(fitted,1)]
        gaps = []
        for r in range(107):
            score = sum(a*row[r] for a,row in zip(residual,rows))
            if r in active:
                assert score.contains(0) and score.rad() < arb("1e-60")
            else:
                gaps.append(-score)
        margin = arb(".00235" if sign == -1 else ".0000235")
        assert min(gaps) > margin
        bound = sum(v*v for v in fitted)
        assert arb(lo) < bound < arb(hi), (bound, sign)
        assert (sign*fitted[-1]-bound).contains(0)
        assert (fitted[-1]**2/sum(v*v for v in fitted)-bound).contains(0)
        print("末阶符号、中心均匀段、最大份额、非活跃半空间最小余量",
              sign,active,bound.str(19),min(gaps).str(15))
        if sign == -1:
            shares = [v/sum(c) for v in c]
            assert shares[-1]/213 > arb(".003")
            print("达到全局最大值的三段份额", *(v.str(17) for v in shares))
    # 复核指定十四项权的实际末圈份额。
    fixed = run_path(str(root/"ch39-height-independent-flavor-vacuum.py"))
    reference = run_path(str(root/"ch39-mincut-reference-ckm-action.py"))
    flavor = run_path(str(root/"ch39-chiral-neutral-flavor-phase.py"))
    winding = run_path(str(root/"ch39-mincut-global-winding.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 224
    x = arb("-.01046")
    y = x-reference["reference_angles"](flavor)[0]
    weights = fixed["minimum_order_weights"](winding,winding["terms"](shell),x,y)
    assert len(weights) == 14 and weights[105] > arb(".072")
    print("指定十四项权的末圈份额",weights[105].str(16))


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
