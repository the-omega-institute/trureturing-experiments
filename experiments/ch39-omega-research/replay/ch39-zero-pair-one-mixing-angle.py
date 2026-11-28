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
    a,b,c,d = direct
    e,f,g,h = inverse
    # (e+gy)(b+dy) - (a+cy)(f+hy) 取零，再用正迹还原 x。
    p0,p1,p2 = e*b-a*f, e*d+g*b-a*h-c*f, g*d-c*h
    discriminant = p1**2-4*p2*p0
    assert discriminant > 0 and abs(p2) > 0
    for sign in (-1,1):
        y = (-p1+sign*discriminant.sqrt())/(2*p2)
        if 0 < y < 1 and abs(b+d*y) > 0 and 0 < solve_x(y) < 1:
            return solve_x(y), y
    raise AssertionError("未能构造六个质量中心的实混合矩阵")


def verify():
    replay = Path(__file__).resolve().parent
    model = run_path(str(replay/"ch39-two-zero-quark-window.py"))
    at_zero = run_path(str(replay/"ch39-quadratic-gram-source.py"))["at_zero"]
    mass_box = [v+arb(0,w.upper()) for v,w in zip(model["ref"][:6],model["tol"][:6],strict=True)]
    u,d = [mass_weights(mass_box[i:i+3]) for i in (0,3)]
    angle = model["ref"][8]+arb(0,model["tol"][8].upper())
    assert (d[2]/u[0])*angle.sin()**2 > arb("3090.46")

    sites = range(-64,65)
    up,down = [{i:at_zero(model,q,n+i) for i in sites}
               for q,n in ((-4,2345666),(2,964890))]
    def normalize(a):
        return a/abs(a.det())**(arb(1)/3)
    direct_down = {i:normalize(a) for i,a in down.items()}
    inverse_up = {i:normalize(a).inv() for i,a in up.items()}
    direct_up = {i:normalize(a) for i,a in up.items()}
    inverse_down = {i:normalize(a).inv() for i,a in down.items()}
    # 质量从轻到重。目标矩阵 D_u V D_d^{-1} 的两因子的奇异值均按降序排列。
    a = [x.sqrt() for x in mass_weights(model["ref"][:3])[::-1]]
    b = [1/x.sqrt() for x in mass_weights(model["ref"][3:6])]
    inverse_a,inverse_b = [[1/x for x in side[::-1]] for side in (a,b)]
    minima = [sum((first[k]*second[2-k])**2 for k in range(3))
              for first,second in ((a,b),(inverse_a,inverse_b))]
    hits = []
    winner = (0,21)
    for i,j in product(sites,repeat=2):
        mat = direct_up[i]*inverse_down[j]
        inverse = direct_down[j]*inverse_up[i]
        tr,inv_tr = [sum(abs(z)**2 for z in m.entries()) for m in (mat,inverse)]
        if (i,j)==winner:
            assert arb("2067.64") < inv_tr < arb("2067.65")
        elif (i,j)!=(0,0):
            assert inv_tr < arb("2067.64")
        if tr < minima[0] or inv_tr < minima[1]:
            continue
        singular = model["singular"](mat)[::-1]
        constraints = ((a,b,singular),(inverse_a,inverse_b,[1/x for x in singular[::-1]]))
        if any(c[k+l]>first[k]*second[l] for first,second,c in constraints
               for k in range(3) for l in range(3-k)):
            continue
        _,y = six_mass_control(model,(tr,inv_tr))
        if (i,j)==(0,-2):
            assert arb(".0000002328") < y < arb(".0000002330")
        hits.append((i,j))
    assert len(hits)==93 and (0,-2) in hits and (0,0) in hits
    print("六个质量中心可达的编号对数",len(hits),"其余 16640 对中逆迹的最大值位于",winner)


if __name__ == "__main__":
    verify()
