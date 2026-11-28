"""实际中央的 10^-12 相位盒，首次连续回返与任意其它盒内正零点的质量差。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


def verify(replay):
    model = run_path(str(replay/"ch39-two-zero-quark-window.py"))
    inverse = run_path(str(replay/"ch39-frozen-three-mass-local-inverse.py"))
    ctx.prec = 384
    a, b = [model["theta"](d)-model["theta"](244) for d in (724,804)]
    ratio = b/a
    eps = arb("1e-12")
    tau = 2*arb.pi()
    threshold = (1+ratio)*eps/tau
    # 格基写作 (724 弦绕数, 804 弦绕数)，两组整数的有向面积为一。
    k0, j0 = 304580939515, 345730428137
    k1, j1 = 1877557160857, 2131218854694
    assert abs(k0*j1-k1*j0) == 1
    e0, e1 = k0*ratio-j0, k1*ratio-j1
    assert e0 > threshold and -threshold < e1 < 0
    # 任意 0 < k < k1 的 (k,j)=u(k0,j0)+v(k1,j1) 必有 u 非零，
    # 且 uv<=0。因此 |k*ratio-j|=|u|e0+|v|(-e1)>=e0。
    for q,n in inverse["CENTERS"].items():
        speed_a, speed_b = abs(q)*a, abs(q)*b
        first = max((tau*k1-eps)/speed_a, (tau*j1-eps)/speed_b)
        last = min((tau*k1+eps)/speed_a, (tau*j1+eps)/speed_b)
        assert first < last
        assert first > (7156521944444 if q == -4 else 14313043888889)
        center = acb.zeta_zero(n).imag
        assert center < first
        # 首次连续回返的整个高度区段在临界线上也能直接排除零点。
        heights = center + first.union(last)
        with ctx.workprec(80):
            modulus = abs(acb(arb(1)/2, heights).zeta())
            assert modulus > arb("1.36" if q == -4 else "0.964")
        # 零绕数时只能落入中央的极窄高度区间；左右相邻编号均在其外。
        assert acb.zeta_zero(n-1).imag < center-eps/speed_a
        assert acb.zeta_zero(n+1).imag > center+eps/speed_a
        offset = model["shell"]["soft"]["point"]["study"]["offset"]
        tail_loss = (arb(23)/25)*arb(2)**-20 / ((center+first)**2+offset)
        assert tail_loss < arb("1.72e-32" if q == -4 else "4.29e-33")
        actual, evaluate = inverse["mass_map"](model,q)
        _, J = evaluate(arb(0,eps.upper()),arb(0,eps.upper()),arb(0).union(actual))
        inverse_jacobian = J.inv()
        bound = (actual-tail_loss)/sum(abs(inverse_jacobian[2,j]) for j in range(3)).upper()
        assert bound > arb("1.99100e-17" if q == -4 else "9.28911e-17")
        print(q,"首个非零连续回返的高度差",first.str(34),"到",last.str(34),
              "宽",(last-first).str(13),"该段 ζ 模",modulus.str(7),"其它盒内实际零点亏额上界",tail_loss.str(13),
              "三质量最大归一差的下界",bound.str(14))


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
