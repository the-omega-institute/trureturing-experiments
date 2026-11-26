"""公共中心商的五态味环、欧氏阈值与热错配界。"""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    source = run_path(str(Path(__file__).with_name("ch39-fermion-arithmetic-prior.py")))
    angle, energy = source["angle"], source["energy"]
    # 统一归一的超荷与两个非阿贝尔分量，单位为一份 Dirac 多重态。
    bidoublet = (F(4,3)*F(3,5)*4*F(1,2)**2, F(4,3)*2*F(1,2), F(0))
    sextet = (F(4,3)*F(3,5)*6*F(1,3)**2, F(0), F(4,3)*2*F(1,2))
    assert tuple(b+c for b,c in zip(bidoublet,sextet)) == (F(4,3),)*3
    assert max(bidoublet)-min(bidoublet) == F(4,3)
    assert F(-16,3)+6*F(2,3) == F(-4,3)
    assert F(-16,3)+4*F(2,3) == F(-8,3)

    # 每个表示的 Dynkin 指数与二次 Casimir 均按 U(1)、SU(2)、SU(3) 排列。
    indices = [(F(3,5), F(1), F(0)), (F(2,5), F(0), F(1))]
    casimirs = [(F(3,20), F(3,4), F(0)), (F(1,15), F(0), F(4,3))]
    adjoint = (0,2,3)
    matrix = [[4*sum(index[i]*casimir[j] for index,casimir in zip(indices,casimirs))
               + (F(20,3)*adjoint[i] if i==j else 0)
               for j in range(3)] for i in range(3)]
    monitor = (5,-3,-2)
    assert all(sum(m*t for m,t in zip(monitor,index)) == 0 for index in indices)
    assert [sum(m*row[j] for m,row in zip(monitor,matrix)) for j in range(3)] == [0,-40,-40]
    ctx.prec = 180
    power = 23226
    target = abs(power*angle(11431)-2*arb.pi()/3)
    second = abs(power*angle(1144)-14*arb.pi()/3)
    margin = arb("2e-6")
    def threshold_bound(delta):
        return -2*(delta/2).cos().log()/(3*arb.pi())
    assert threshold_bound(target) < arb("1.811e-11")
    assert threshold_bound(target+margin) < arb("2.099e-11")
    # 双双态使费米子取向能加倍，六重态的质量矩阵不含取向。
    effective_exponent = arb("2e15")
    kappa = arb("1e23")
    single = margin/(power-1)
    tail = 22*arb.pi()*kappa.sqrt()*(arb(1)/2-2*kappa*(single/2).sin()**2).exp()
    gap = energy(second-margin)-energy(target+margin)
    cost = energy(target+margin)-energy(arb(0))
    weight_ratio = (1+arb.const_euler()/2-(4*arb.pi()).log()/2)*(arb(1)/4+acb.zeta_zero(11431).imag**2)
    # 在坏角事件上纳入全部编号和取向，在好角事件上只纳入非目标。
    failure = (2*weight_ratio/(1-tail)*
               ((-effective_exponent*gap).exp()+tail*(effective_exponent*cost).exp()))
    assert failure < arb("6.03e-107")
    print("双双态与六重态的一圈系数",bidoublet,sextet,
          "中心和好角窗的三支差分阈值上界", threshold_bound(target).str(12),
          threshold_bound(target+margin).str(12),"联合错配上界", failure.str(12))


if __name__ == "__main__":
    verify()
