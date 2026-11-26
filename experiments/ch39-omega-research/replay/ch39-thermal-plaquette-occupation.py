"""有限零点与计次集合的 Gibbs 权重：同时积分连续尺度并比较全部长度。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root=Path(__file__).parent
chain=run_path(str(root/"ch39-local-winding-chain.py"))
joint=chain["joint"]
q=arb


def verify():
    ctx.prec=128
    count=23226
    theta=joint["winding"]["angle"](11431)
    alpha=2*q.pi()/3
    delta=alpha-count*theta
    energy=lambda a: 2*(1-a.cos())
    slope=3/(1+q(3)/8)**2
    early_drop=4*(theta/2).sin()*min((delta+theta/2).sin().lower(),alpha.sin().lower())
    early=slope*early_drop/count-q("2e-16")
    late=slope*(energy(theta-delta)/(count+21000)-energy(delta)/count)
    assert theta-delta<2*alpha-(count+21000)*theta
    tail=21000*q("1e-17")-3*energy(delta)/count
    temperature=q("2e-17")
    log_sum=21001
    assert q("1e-17")/temperature>q(".49")
    # 三段分别给几何级数、有限求和、几何级数的上界。
    smallest=min((early/temperature).lower(),(late/temperature).lower(),(tail/temperature).lower())
    p_exponent=smallest-q(log_sum+10).log()
    print("other counts: negative log upper bound",p_exponent.str(12),flush=True)
    assert early/temperature > 1 and p_exponent>4000

    time=q("26.81676714125598")
    width=q("1e-8")
    target=joint["flow"]["gauge"]["shortest"]["model"]["target"]
    center=target(11431,time)
    allowance=q("1.005e-7")+(q(".09")+2*center*(2*width).exp())*width
    spread=q(".1435")*q(".00005")-2*q("1.005e-7")
    difference=min((target(11430,time+q(".00005"))-center).lower(),
                   (center-target(11432,time-q(".00005"))).lower())-q("1.005e-7")-q(".09")*q(".00005")
    trial=3*energy(delta)/count+count*q("2e-16")
    for leads in (1,2,3):
        m,d=[q(v) for v in joint["leads"]["terminal_modes"](leads)]
        competitor=min((3*m*difference**2).lower(),(d*spread**2/2).lower())
        gap=competitor-3*d*allowance**2-trial
        # 每个竞争编号的尺度区间总长 10；试探扇区至少积分宽 2w。
        # 1/(1-exp(-.49)) < 3。
        assert 1/(1-(-q(".49")).exp())<3
        exponent=gap/temperature-(15/width).log()
        print("lead count",leads,"other label, negative log upper bound",exponent.str(12),flush=True)
        assert exponent>200000
    height=acb.zeta_zero(11431).imag
    norm=1+q.const_euler()/2-(4*q.pi()).log()/2
    assert (norm*(q(1)/4+height**2)).log()<16


if __name__=="__main__":
    verify()
