"""四个活跃钟码和两组对偶不等式确定混合味读数的全局最低权重。"""
from math import comb
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, arb_poly, ctx

previous = run_path(str(Path(__file__).with_name("ch39-mixed-flavor-pivot-weight.py")))
support = (11431,11831,11032,11432)


def determinant(m):
    if len(m)==1:
        return m[0][0]
    return sum((-1)**k*h*determinant([v[:k]+v[k+1:] for v in m[1:]])
               for k,h in enumerate(m[0]))


def positive_on_unit(p):
    n = p.degree()
    assert n >= 0
    assert all(sum(p[j]*comb(k,j)/comb(n,j) for j in range(k+1))>0 for k in range(n+1))


def verify(source=None):
    ctx.prec = 224
    points = previous["points"] if source is None else source
    rows = {i:tuple(map(arb,z)) for i,z in points()}
    assert set(rows) == set(range(10815,12093))
    S,T = [arb((previous[key]-previous["half"])**2) for key in ("s12","s23")]
    lo, hi = [arb((previous["s13"]+k*previous["half"])**2) for k in (-1,1)]
    tau = arb("1.141").tan()
    U = arb_poly([lo,hi-lo])
    def columns(row):
        a,b,u,d,j = row
        dlow = ((1-U)*d+(-1+(1+U)*T)*a+(-U+(1+U)*S)*b-(1-U**2)*S*T)*arb('0.5')
        return [arb_poly([1]),arb_poly([a-S*(1-u)]),arb_poly([b-T*(1-u)]),j-tau*dlow]
    m = [columns(rows[i]) for i in support]
    det = determinant(m)
    sign = 1 if det(0)>0 else -1
    det *= sign
    positive_on_unit(det)
    for mu, direction in ((0,-1),(55000,1)):
        rhs = [arb_poly([int(i==11431)-mu*rows[i][2]]) for i in support]
        num = [sign*determinant([v[:k]+[h]+v[k+1:] for v,h in zip(m,rhs)]) for k in range(4)]
        for p in num[1:]:
            positive_on_unit(p)
        derivative = mu*(hi-lo)*det**2 + num[0].derivative()*det - num[0]*det.derivative()
        positive_on_unit(direction*derivative)
        for n,z in rows.items():
            if n in support:
                continue
            slack = (int(n==11431)-mu*z[2])*det - sum(p*x for p,x in zip(num,columns(z)))
            positive_on_unit(slack)
    # 对四个支撑点同时指定总重、两条角下边及均值 U，得到仿射权重。
    k = arb_mat([[1]*4,
                 [rows[i][0]-S*(1-rows[i][2]) for i in support],
                 [rows[i][1]-T*(1-rows[i][2]) for i in support],
                 [rows[i][2] for i in support]])
    intercept = k.solve(arb_mat([[1],[0],[0],[lo]]))
    slope = k.solve(arb_mat([[0],[0],[0],[hi-lo]]))
    weights = [arb_poly([intercept[i,0],slope[i,0]]) for i in range(4)]
    d,j = [sum(p*rows[i][c] for p,i in zip(weights,support)) for c in (3,4)]
    phase = 2*j-tau*(1-U)*(d-S*(1-T)-(1-S)*T*U)
    assert phase(0)*phase(1)<0
    positive_on_unit(phase.derivative() if phase(1)>0 else -phase.derivative())
    left,right = arb(0),arb(1)
    left_sign = phase(left)>0
    for _ in range(90):
        middle = (left+right)/2
        val = phase(middle)
        assert val>0 or val<0
        if (val>0)==left_sign:
            left=middle
        else:
            right=middle
    root = left.union(right)
    answer = [p(root) for p in weights]
    assert all(w>0 for w in answer) and j(root)>0
    if source is None:
        assert arb(".4675306")<answer[0] and answer[0]<arb(".4675308")
    print("Unique minimizer on closed flavor boxes",*(w.str(15) for w in answer),
          "averaged |Vub|²",U(root).str(15))


if __name__ == "__main__":
    verify()
