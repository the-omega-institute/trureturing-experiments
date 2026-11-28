"""第120与163枚正零点共用八站电流时的近最优振幅界。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


def verify(replay):
    laplacian = run_path(str(replay/"ch39-three-ray-chiral-boundary.py"))["laplacian"]
    ctx.prec = 256
    sites = (0,2,3,5,6,9,11,13)
    greens = [laplacian(5, acb.zeta_zero(n).imag).inv() for n in (120,163)]
    r, t = arb(".08216966"), arb(".54905026")
    upper = arb(".00003295347")
    for green in greens:
        amplitude = ((1-r)*(t*green[6,3]+(1-t)*green[6,5])
                     + r*(t*green[0,3]+(1-t)*green[0,5]))
        assert abs(amplitude) < upper
        print("两端共同见证的振幅模长", abs(amplitude).str(16))

    for i in sites:
        for j in sites:
            special = {i,j} in ({6,3},{6,5})
            assert greens[1][i,j].real > (arb("-.000532") if special else arb(".006"))
    k = arb(".909")
    c = (arb(".006")-upper)/arb(".006532")
    assert c-(1-c)**2/(2*c) > k

    def vertex(i, j):
        return tuple(k*g[6,j]+(1-k)*g[i,j] for g in greens)

    cases = (
        ("0", ".54912", (-336310,359998,368088,-349165), ".000032950060"),
        (".54912", "1", (-336684,353531,368100,-355585), ".000032950141"),
    )
    for lo, hi, weights, threshold in cases:
        x,y,u,v = map(arb, weights)
        denominator = (x*x+y*y).sqrt()+(u*u+v*v).sqrt()
        bounds = []
        for split in (arb(lo), arb(hi)):
            for i in sites:
                inner = tuple(split*a+(1-split)*b for a,b in zip(vertex(i,3),vertex(i,5)))
                pairs = [inner]
                for j in sites:
                    if j not in (3,5):
                        pairs.append(tuple(k*a+(1-k)*b for a,b in zip(inner,vertex(i,j))))
                for a,b in pairs:
                    projection = (x*a.real+y*a.imag+u*b.real+v*b.imag)/denominator
                    assert projection > arb(threshold)
                    bounds.append(projection.lower())
        assert arb(threshold) > arb(".0000329500")
        print("接收站3相对权重的区间与归一对偶下界", lo, hi, min(bounds).str(15))


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
