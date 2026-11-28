"""在两编号的 129×129 邻域检验仅用六质量及最小混合角能否辨识中央零点。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb


def mass_weights(masses):
    return [(v**2/(masses[(i+1)%3]*masses[(i+2)%3]))**(arb(2)/3)
            for i,v in enumerate(masses)]


def six_mass_control(model, traces):
    # θ13=0 时 x=sin² θ12、y=sin² θ23 的四个双线性系数。
    up,down = [mass_weights(model["ref"][i:i+3]) for i in (0,3)]
    def coefficients(s,target):
        w = [[(u/d)**s for d in down] for u in up]
        return (w[0][0]+w[1][1]+w[2][2]-target,
                -w[0][0]+w[0][1]+w[1][0]-w[1][1],
                -w[1][1]+w[1][2]+w[2][1]-w[2][2],
                -w[1][0]+w[1][1]+w[2][0]-w[2][1])
    direct,inverse = [coefficients(s,t) for s,t in zip((1,-1),traces,strict=True)]
    def solve_x(y):
        a,b,c,d = direct
        return -(a+c*y)/(b+d*y)
    def residual(y):
        a,b,c,d = inverse
        return a+c*y+(b+d*y)*solve_x(y)
    left,right = arb(".0000002328"),arb(".0000002330")
    assert residual(left)*residual(right) < 0
    y = left.union(right)
    denominator = direct[1]+direct[3]*y
    assert not denominator.contains(0)
    assert 0 < solve_x(y) < 1 and 0 < y < 1


def verify():
    replay = Path(__file__).resolve().parent
    model = run_path(str(replay/"ch39-two-zero-quark-window.py"))
    at_zero = run_path(str(replay/"ch39-quadratic-gram-source.py"))["at_zero"]
    normalized_traces = run_path(str(replay/"ch39-two-zero-sharp-quark-fit.py"))["normalized_traces"]
    mass_box = [v+arb(0,w.upper()) for v,w in zip(model["ref"][:6],model["tol"][:6],strict=True)]
    u,d = [mass_weights(mass_box[i:i+3]) for i in (0,3)]
    angle = model["ref"][8]+arb(0,model["tol"][8].upper())
    assert (d[2]/u[0])*angle.sin()**2 > arb("3090.46")

    sites = range(-64,65)
    up,down = [{i:at_zero(model,q,n+i) for i in sites}
               for q,n in ((-4,2345666),(2,964890))]
    six_mass_control(model, normalized_traces(up[0], down[-2]))
    def normalize(a):
        return a/abs(a.det())**(arb(1)/3)
    direct_down = {i:normalize(a) for i,a in down.items()}
    inverse_up = {i:normalize(a).inv() for i,a in up.items()}
    winner = (0,21)
    for i,j in product(sites,repeat=2):
        if (i,j)==(0,0):
            continue
        matrix = direct_down[j]*inverse_up[i]
        value = sum(abs(z)**2 for z in matrix.entries())
        if (i,j)==winner:
            assert arb("2067.64") < value < arb("2067.65")
        else:
            assert value < arb("2067.64")
    print("其余 16640 对中逆迹的最大值位于编号偏移", winner, "并小于 2067.65")


if __name__ == "__main__":
    verify()
