"""按三路径核心的四个独立顶点相位投影量子热迹。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

root = Path(__file__).parent
q = arb


def verify():
    ctx.prec = 128
    quantum = run_path(str(root/"ch39-quantum-endpoint-chain.py"))
    volumes = run_path(str(root/"ch39-classical-endpoint-ensemble.py"))["volume_bounds"]()
    p = 23226
    heat, hmin = q("2e-17"),q("5e-18")
    theta = run_path(str(root/"ch39-first-arithmetic-plaquette-winding.py"))["angle"](11431)
    alpha = 2*q.pi()/3
    e = 2*(1-(p*theta-alpha).cos())/p
    reverse = 2*(1-(p*theta+alpha).cos())/p
    slope = 3*(reverse-e)/(2*(1+e/8)*(1+reverse/8))
    for kmax in (2,10000):
        thermal_span = volumes+quantum["packet_cost"](q(kmax),q(1))
        # Jensen 的平均能量上界。F(H)/H-F(2H)/(2H) 为逐频率 log(2cosh(u/2)) 之和。
        energy = 3*e+2*heat*(10*q(2).log()+thermal_span+6*(p-1)*(2*q(2).cosh()).log())
        momentum_sq = 4*energy/hmin**2
        radius_sq = 9+12*p+(q(p)**2+q(p)/2+q(1)/8)*energy
        charge = (momentum_sq*radius_sq).sqrt()
        assert charge > 1
        # 角四维立方体内部扭曲热迹 >= Z/2，余下扭角的实热核迹非负。
        loss = q(2).log()+4*(8*q.pi()*charge).log()
        assert loss < 154
        power = 4000-thermal_span-loss-q(2).log()
        chi_deficit = 2*heat*(thermal_span+loss+q(16).log())/slope
        print("闭合系数上端",kmax,"平均总场能上界",energy.str(9),
              "投影损失负对数上界",loss.str(9),
              "错配两项中较小指数的下界",power.str(9),"取向缺额上界",chi_deficit.str(10))
        if kmax==2:
            assert power>3729 and chi_deficit < q("5.7e-11")
        else:
            assert power>1423 and chi_deficit < q("5.4e-10")


if __name__ == "__main__":
    verify()
