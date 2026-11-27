"""检验固定零点盒的六质量是否允许共同源协方差有重本征值。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def verify():
    orbit = run_path(str(Path(__file__).with_name("ch39-two-zero-orbit-condition.py")))
    model,up,down = orbit["calibration"]()
    adj = model["adj"]
    g_up,g_down = [adj(a)*a/abs(a.det())**(arb(2)/3) for a in (up,down)]
    (g,u),(h,d) = [model["eig"](a) for a in (g_up,g_down)]
    smallest_g = 1/(g[1]*g[2])
    smallest_h = 1/(h[1]*h[2])
    assert smallest_g < arb(".000157") and smallest_h < arb(".000557")
    assert g[1] < arb("12.728") and g[2] < 501 and h[2] > 477
    assert g_up.trace().real < 514 and 1/smallest_g+1/g[1]+1/g[2] < 6380

    identity = acb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    difference = 4*(g_down.trace().real*identity-g_down)-(g_up.trace().real*identity-g_up)
    for n in (1,2,3):
        minor = acb_mat([[difference[i,j] for j in range(n)] for i in range(n)])
        assert minor.det().real > 0

    masses = [x+arb(0,r.upper()) for x,r in zip(model["ref"][:6],model["tol"][:6],strict=True)]
    def traces(values):
        scale = (values[0]*values[1]*values[2])**(arb(2)/3)
        return sum(x*x for x in values)/scale,scale*sum(1/(x*x) for x in values)
    fu,iu = traces(masses[:3])
    fd,idown = traces(masses[3:])
    assert fu > 110000 and 158000 < iu < 185000 and fd < 1485 and idown < 794
    assert fu/fd > 74

    overlap = adj(u)*d
    for i,j,limit in ((0,0,".834"),(1,0,".554"),(2,0,".0125"),(2,1,".0161")):
        assert abs(overlap[i,j]) < arb(limit)
    assert (1-abs(overlap[0,0])**2).sqrt() < arb(".555")

    lower,upper = arb(".0083"),arb(".041")
    assert arb(6380)/158000 < upper
    assert 794*smallest_h < arb(".443") and 185000*smallest_g < arb("29.1")
    assert (arb(".834")+arb(".555")*(arb(".443")*lower).sqrt())**2 < 1-arb("29.1")*lower
    assert (arb(1485)/477).sqrt() < arb("1.77")
    transverse = (arb(".443")*upper).sqrt()
    assert transverse < arb(".135")
    assert arb(".0125")+arb(".0161")*transverse < arb(".0147")
    assert smallest_g+g[1]*(arb(".554") + transverse)**2 < arb("6.1")
    # 对 lambda>lower，右侧除以 lambda² 后的两项均递减。
    assert arb("6.1")/lower**2 + 501*(arb(".0147")/lower+arb("1.77"))**2 < 95000
    assert 514*upper+95000 < 95100
    print("重根 λ ≥ 1 的归一化迹比上界",4)
    print("重根 λ < 1 时的强制区间、上型正迹上界、数据下界",lower.str(5),upper.str(5),95100,110000)


if __name__ == "__main__":
    verify()
