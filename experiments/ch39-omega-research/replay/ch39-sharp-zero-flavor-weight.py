"""直接输入真实零点高度和解析球壳响应时的混合味最低权重。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent
sharp = run_path(str(root/"ch39-sharp-flavor-pivot-weight.py"))


def observables():
    # 在设置精度之后装入带预计算切平面基的原模型。
    flavor = run_path(str(root/"ch39-fixed-local-flavor.py"))
    coeff = [arb(n)/2048 for n in flavor["codes"]]
    for n,zero in enumerate(acb.zeta_zeros(10815,1278),10815):
        w = flavor["mixing"](flavor["native"]["plane_response"](zero.imag),coeff)
        j = (w[0,0]*w[1,0].conjugate()*w[1,1]*w[0,1].conjugate()).imag
        yield n, [abs(w[i,k])**2 for i,k in ((0,1),(1,2),(0,2),(1,0))]+[j]


if __name__ == "__main__":
    sharp["verify"](observables)
