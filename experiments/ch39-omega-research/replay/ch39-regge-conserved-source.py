"""从 Regge 平直胞的长波 Hessian 复核守恒源极点和洛伦兹零锥的秩。"""
from pathlib import Path
from runpy import run_path

import sympy as sp


def verify():
    root = Path(__file__).parent
    fp = run_path(str(root/"ch39-regge-fierz-pauli.py"))
    fp["base"]["verify"]()
    fp["verify"]()
    metric = sp.eye(4)-sp.ones(4)/5
    inverse = metric.inv()
    k = sp.Matrix(sp.symbols("k0:4"))
    raised = inverse*k
    momentum_squared = (k.T*raised)[0]
    basis = []
    for i,j in fp["METRIC"]:
        h = sp.zeros(4)
        h[i,j] = h[j,i] = 1
        basis.append(h)
    trace = sp.Matrix([(inverse*h).trace() for h in basis])
    gram = sp.Matrix([[(inverse*h*inverse*v).trace() for v in basis] for h in basis])
    gauge = sp.Matrix.hstack(*(h*raised-k*tr/2 for h,tr in zip(basis,trace)))
    dual = gram-trace*trace.T/2
    fixed = fp["continuum"](list(k))+2*gauge.T*inverse*gauge
    assert (fixed-momentum_squared*dual).applyfunc(sp.expand) == sp.zeros(10)
    assert dual*gram.inv()*dual == gram
    xi = sp.Matrix(sp.symbols("u0:4"))
    shift = k*xi.T+xi*k.T
    assert (shift*raised-k*(inverse*shift).trace()/2-momentum_squared*xi).applyfunc(sp.expand) == sp.zeros(4,1)

    # 非零零向量沿第三空间轴。守恒实对称源的六个独立分量。
    a,b,c,p,m,n = sp.symbols("a b c p m n",real=True)
    source = sp.Matrix([[a,b,c,a],[b,p,m,b],[c,m,n,c],[a,b,c,a]])
    lorentz = sp.diag(-1,1,1,1)
    assert source*sp.Matrix([-1,0,0,1]) == sp.zeros(4,1)
    numerator = (lorentz*source*lorentz*source).trace()-(lorentz*source).trace()**2/2
    assert sp.expand(numerator-(p-n)**2/2-2*m*m) == 0
    assert sp.hessian(numerator,(a,b,c,p,m,n)).rank() == 2
    print("归一 Hessian 的守恒源极点系数、非零零动量壳上的留数秩、静态单位源系数",
          20,2,10)


if __name__ == "__main__":
    verify()
