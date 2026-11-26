"""在减法点分类能补偿整个中性交换的秩一直接匹配，并检验第二个动量。"""
from pathlib import Path
from runpy import run_path

from flint import arb_mat

light = run_path(str(Path(__file__).with_name("ch39-light-weyl-current-basis.py")))
poles, mat, scheme, q = (light[key] for key in ("poles", "mat", "scheme", "q"))


def axes(e, z):
    c = poles["colour"](-e, poles["unit"])/24
    r = poles["right"](-e, poles["unit"])/4
    s = 9*c+r
    h = z+s
    g = arb_mat([[3*c], [-r], [0]])
    v = arb_mat([[-6*c/h], [2*r/h], [1]])
    alpha = 36*e*c*r/s
    beta = e*z/(s*h)
    assert min(c, r, alpha, beta) > 0
    return h, g, v, alpha, beta


def same_matrix(candidate, reference):
    residual = candidate-reference
    assert all(x.contains(0) and abs(x)<q("1e-24") for x in residual.entries())


def verify():
    time = poles["study"]["nominal"][2]
    targets = [poles["study"]["comp"]["model"]["target"](n, time) for n in (11430, 11432)]
    for i, logs in enumerate(poles["retune"]["log_sets"]):
        (length, z, rho, raw_d, raw_b, _), _, _ = poles["retune"]["predictions"](logs, time)
        e = scheme["anchor"]/length
        b, d, k = raw_b/rho, raw_d/rho, 1/rho
        h = mat["flow"]["running"]["hypercharge_window"](time+scheme["anchor"].log()/2).mid()**2
        side, g, v, alpha, beta = axes(e, z)
        p = k*e/(b+d*e+4/side)
        coefficient = h+k*(alpha-q("1/24"))
        bare = mat["massive_exchange"](e,z,b,d)*k+mat["hypercharge"]*h
        gg = g*g.transpose()
        same_matrix(mat["hypercharge"]*coefficient+gg*(k*beta)+v*v.transpose()*p, bare)
        for target, sign in zip(targets, (1, -1), strict=True):
            kk = 1/(poles["offset"]+target/(15*length))
            hh = coefficient-kk*(alpha-q("1/24"))
            dd = (kk*e/p-4/side-b)/e
            assert 0 < hh < 1 and dd > 0
            update = (k-kk)*beta
            candidate = mat["massive_exchange"](e,z,b,dd)*kk+mat["hypercharge"]*hh
            assert mat["positive"](candidate)
            same_matrix(candidate+gg*update,bare)
            assert sign*update > 0
        g2 = axes(q("5/3")*e,z)[1]
        ratio = g[1,0]/g[0,0]
        ratio2 = g2[1,0]/g2[0,0]
        assert q(".0195276") < abs(ratio2-ratio) < q(".0195279")
        print("case", i, "R/D slopes", ratio.str(13), ratio2.str(13),
              "absolute drift", abs(ratio2-ratio).str(10), flush=True)


if __name__ == "__main__":
    verify()
