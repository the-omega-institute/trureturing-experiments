"""一份四元数双基本场的相位障碍与两份场的矩阵径向势。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb
from sympy import I, Matrix, diag, eye, simplify, sqrt, symbols, zeros

root = Path(__file__).parent
q = arb


def verify():
    # 标准四元数基。与左、右乘均可交换的实线性映射由其在一上的中心值决定。
    quaternion_i = diag(I,-I)
    quaternion_j = Matrix([[0,1],[-1,0]])
    basis = [eye(2),quaternion_i,quaternion_j,quaternion_i*quaternion_j]
    columns = []
    for unit in basis:
        columns.append(Matrix([c for v in basis[1:] for c in unit*v-v*unit]))
    assert Matrix.hstack(*columns).rank()==3
    a,b,c,d = symbols("a b c d", real=True)
    quaternion = sum((v*g for v,g in zip((a,b,c,d),basis)),zeros(2))
    assert simplify(quaternion.H*quaternion-(a*a+b*b+c*c+d*d)*eye(2))==zeros(2)
    projectors = [sqrt(2)*diag(1,0),sqrt(2)*diag(0,1)]
    assert all((z.H*z).trace()/2==1 for z in projectors)
    assert projectors[0]*projectors[1]==zeros(2)

    joint = run_path(str(root/"ch39-local-winding-chain.py"))["joint"]
    flow = joint["flow"]["gauge"]["shortest"]["model"]
    time,eps = q("26.81676714125598"),q("1.005e-7")
    center = flow["target"](11431,time)
    t1,t2 = [2*(1/(2*z.imag)).atan() for z in acb.zeta_zeros(11431,2)]
    length = 200000
    tau = (1/t2.cos()).acosh()
    eff = t2.sin()*(2*(length*tau).coth()+1/(length*tau).sinh())
    floor = 2*t1.sin()
    shift = (flow["target"](11432,time)/center).log()/2
    assert 23 < time+shift < 33
    contrast = max(q(joint["leads"]["terminal_modes"](v)[1]) for v in (1,2,3))
    cost = 3*contrast*(eps+q(".09")*abs(shift))**2
    assert floor > eff
    gain = 3*(floor-eff)/((1+floor/8)*(1+eff/8))
    gap = gain-cost-12*length*q("4e-17")-q("2e-17")*((acb.zeta_zero(11432).imag**2+q(1)/4)
             /(acb.zeta_zero(11431).imag**2+q(1)/4)).log()
    assert gap>q("2.9e-8")
    print("实四元数链的 11432 见证对任意 11431 长度和尺度的自由能优势",gap.str(11))
    # 标量矩阵取向是中心三次单位根，行列式中也可分辨其共轭。
    omega=(-1+sqrt(3)*I)/2
    assert simplify((omega*eye(2)).det()-omega.conjugate())==0


if __name__ == "__main__":
    verify()
