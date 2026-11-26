"""局域汤川桥命中上下型范数比时真空方向的重质量代价。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

study=run_path(str(Path(__file__).with_name("ch39-heavy-higgs-flavor-fraction.py")))
local=study["local"]
adj=study["adj"]


def verify():
    ctx.prec=384
    eye=acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    ratio=(sum(arb(v)**2 for v in (".00000615",".00311",".8616"))
           /sum(arb(v)**2 for v in (".0000135",".000268",".01401"))).sqrt()
    low,high=ratio*arb(".999")/arb("1.001"),ratio*arb("1.001")/arb(".999")
    unit=arb("2.435e18")/724**2
    for number in (22027,28325):
        h=acb.zeta_zero(number).imag
        lap=local["local"]["shell"]["laplacian"]
        neutral=lap(0,h)
        center=sum(neutral[i,i].real for i in range(19))/18
        V={q:local["basis"](q,h,center) for q in (1,-4,2,3)}
        A=study["operator"](V[1],V[-4])
        D=study["operator"](V[1],V[2]).conjugate()
        B=adj(A)*A/low**2-adj(D)*D
        P=eye-V[3]*adj(V[3])
        L=lap(3,h)
        two_step=L*(L-center*eye)
        # H₀ 是完整质量平方算子。给予前三个轻方向零质量，使它们的代价最小。
        mass_squared=P*two_step**2*P
        normalization=arb(10000)
        rayleigh=mass_squared/normalization
        rl,rh=study["certificate"](rayleigh,B,eye,A,D,low,high)
        mh2=(arb("125.25")/unit)**2
        residual=((mass_squared-mh2*eye)/normalization)**2
        vl,vh=study["certificate"](residual,B,eye,A,D,low,high)
        ray_low,ray_high=unit*(normalization*rl).sqrt(),unit*(normalization*rh).sqrt()
        norm_low,norm_high=unit**2*normalization*vl.sqrt(),unit**2*normalization*vh.sqrt()
        print(number,"最少 Rayleigh 质量/GeV 与 Hermitian 质量平方改动/GeV² 的包围",
              *(v.str(8) for v in (ray_low,ray_high,norm_low,norm_high)))
        assert ray_low>arb("1e13") and norm_low>arb("1e27")


if __name__=="__main__":
    verify()
