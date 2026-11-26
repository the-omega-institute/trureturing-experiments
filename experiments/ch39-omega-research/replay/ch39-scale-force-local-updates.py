"""以规范能在两个端点的尺度导数控制可逆局部跃迁。"""
from pathlib import Path
from runpy import run_path

from flint import arb

conditional = run_path(str(Path(__file__).with_name("ch39-conditional-scale-resonance.py")))
q = arb


def verify():
    shifts, reference, mean, pairs = conditional["verify"]()
    joint = conditional["source"]["joint"]
    large = max(abs(v).upper() for v in shifts)
    step_b = max(abs(b-a).upper() for a,b in zip(shifts,shifts[1:]))
    eps, heat = q("1.005e-7"), q("2e-17")
    fine, coarse = q("2e-8"), q(".001")
    reach, force, gate_core = q(".0005"), q("1.5e-6"), q("2e-8")
    curvature_low, curvature_high = 20, 50
    assert large + q(".01") < q(".07") and 50*gate_core < force
    for j in (1,2,3):
        m,d = (q(v) for v in joint["leads"]["terminal_modes"](j))
        # 带外两侧共同端口残差的绝对值下界。
        far = reference*(1-(-q(".02")).exp())-eps-mean.upper()*(large+q(".01"))
        assert far > 0
        assert 6*m*far*mean.lower() > q(".0001")
        assert all(v.lower()*q(".07")>2*eps for v in pairs)
        near_slope = 2*reference*(-2*(large+q(".07"))).exp()+mean.lower()
        contrast_rate = 2*d*sum((2*eps+q(".07")*v.upper())*v.upper() for v in pairs)/3
        assert 6*m*far*near_slope-contrast_rate > q(".0001")
        # 对连续 b 的隐函数斜率作包络，无需求出实际的 1278 个极小点。
        qmax = reference*q(".02").exp()
        residual = eps + mean.upper()*(large+q(".01")) + reference*(q(".02").exp()-1)
        mixed = 6*m*(2*qmax*(mean.upper()+2*qmax)+4*qmax*residual)
        radius = mixed/curvature_low*step_b
        assert radius+force/curvature_low+gate_core < reach
    assert q(".0001") > force
    assert curvature_low*fine < q(".0001")
    assert curvature_low*fine**2/(2*heat) > 190
    inside_excess = q(curvature_high)/2*(force/curvature_low)**2
    assert (q(20000)*heat-inside_excess)/heat > 12000
    max_steps, warmup = 300000000, 1000000
    def lower_tail(mean_count, required):
        assert mean_count > required
        return (mean_count-required)**2/(2*mean_count)
    funds = [lower_tail(q(200000)/8, 10001), lower_tail(q(800000)/8,50001),
             lower_tail(q(warmup)/8,24385),
             lower_tail(q(max_steps-warmup)*gate_core/(8*reach),661)]
    risk = (sum((-v).exp() for v in funds)
            +max_steps*((-q(64)).exp()+(-q(190)).exp()+(-q(12000)).exp()))
    assert risk < q("5e-20")
    print("局部提议步数",max_steps,"末态失败概率上界",risk.upper().str(12))


if __name__ == "__main__":
    verify()
