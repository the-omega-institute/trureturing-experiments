"""一个跨部门的原始一阶迹差何时把两编号同时限制到有限范围？"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    replay = Path(__file__).resolve().parent
    shell = run_path(str(replay/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 320
    shift = shell["soft"]["point"]["study"]["offset"]
    scale = 6*arb(23)/25*arb(2)**-20
    def missing(n):
        h = acb.zeta_zero(n).imag
        return scale/(h*h+shift)

    def inverse(y):
        h = (scale/y-shift).sqrt()
        count = h.zeta_nzeros().unique_fmpz()
        assert h > 0 and count is not None
        n = int(count)
        assert acb.zeta_zero(n).imag < h < acb.zeta_zero(n+1).imag
        return n

    target = missing(964890)-missing(2345666)
    k = inverse(target)
    assert k == 1087795
    above = missing(k)-target
    below = target-missing(k+1)
    assert 0 < below < above < target
    # 较近的聚点在目标之下，所有有限上型编号只会继续降低该行的差值。
    half = below/2
    assert arb(2)**-78 < half < arb(2)**-77
    width = 2*arb(2)**-78
    upper_limit = inverse(above-width)
    assert upper_limit == 2367696686
    assert abs(missing(k)-missing(upper_limit)-target) < width
    assert missing(k)-missing(upper_limit+1)-target > width
    print("差值、最近聚点所在下型编号、有限性临界闭半宽",
          target.str(17),k+1,half.str(17))
    print("在绝对半宽 2^-78 下可同时达到的最大上下型编号",upper_limit,k)


if __name__ == "__main__":
    verify()
