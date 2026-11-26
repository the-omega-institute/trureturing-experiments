"""一般二阶矩阵中介的非标量分量在竞争盒内严格增能。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

scalar = run_path(str(Path(__file__).with_name("ch39-free-complex-messengers.py")))
path_energy = run_path(str(Path(__file__).with_name("ch39-matrix-messenger-path.py")))


def verify(energy_slack="0"):
    path_energy["verify"]()
    scalar["verify"](energy_slack)
    ctx.prec = 160
    r = arb("4e-5")
    R = arb("3e-18")*(arb(25000)/42)**4
    u = arb("5e-18")*R
    C = 4/(8*arb.pi()**2)
    # Frobenius 归一模长至多二，三个跳块的最大算子范数小于三。
    tau = 6*r+9*r**2
    remainder = 3*C*tau**3/(6*(1-tau))
    dmin, dmax = arb(4), arb("4.8")
    p = 23226
    theta = scalar["free"]["absolute"]["angle"](11431)
    delta = 2*arb.pi()/3-p*theta
    end_cost = R*12*(delta/2).sin()**2/p
    # 试探取三个单位模的等边相位。一圈能上下各用一次一致余项。
    budget = (2*remainder + 4*(C*r**2)**2/(R*dmin)
              + end_cost + u*3*arb.pi()*(p-1)/2 + arb(energy_slack))
    eps = (budget/R).sqrt()
    assert eps < arb(".001381")
    eps = arb(".001381")
    # 扣除半份原路径能以吸收一圈的标量质量主项。
    sing = 2*eps/dmin.sqrt()
    assert sing < arb(".002")
    lower, upper = 1-sing.upper(), 1+sing.upper()
    # 将 B0 的极分解酉因子吸收进基底。b0 为实正的平均奇异值。
    # Gram 行列式控制 Bj 在 B0 的正交补上的距离。
    tracer0 = sing.upper()
    tracers = (eps.upper()/lower + (upper/lower)*tracer0)
    assert tracers < arb(".0041")
    min_b = (lower**2-tracers**2).sqrt()
    assert min_b > arb(".997")
    # 每条链端的标量投影残差除以 q 至多为 eps + tracers²。
    chain_loss = 2*(eps+tracers**2)
    shape_loss = dmax*(1/min_b-1) + 2*(1-min_b**2) + chain_loss
    # 支撑 d||E0||² 与两块含 B0 的楔积的显式强制性常数。
    gram_factor = min_b**2/2
    anchor_factor = dmin/(1+4*upper**2/min_b**2)
    positive = min(gram_factor.lower(), anchor_factor.lower())
    # 标量投影的一圈梯度在无迹方向为零。只需二阶负部。
    first = 2*arb(2).sqrt()*r*(1+3*r)
    second = 4*r*r
    hess_error = 3*C*((-2*(1-3*r).log())*first**2 + tau**2/(2*(1-tau))*second)
    margin = positive-shape_loss-hess_error/(2*R)
    assert margin > arb(".44")
    print("Tree budget sqrt(T/R)",eps.upper().str(9),"traceless norm",tracers.upper().str(9),
          "energy gain above scalar projection per R sum norm²",margin.lower().str(9))


if __name__ == "__main__":
    verify()
