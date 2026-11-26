"""两种原生钟弦标定三代边界的质量单位，再从两圈电弱交点反求强耦合。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

ROOT = Path(__file__).parent
two = run_path(str(ROOT / "ch39-chiral-center-two-loop.py"))
sm, chiral = two["sm"], two["chiral"]


def gauge_segments(L,possible,steps):
    b,B = sm["BETA"],sm["B"]
    segments = []
    last = arb(0)
    for i,(l,w,extra,_,_) in enumerate(possible):
        end = l.min(L)
        dt = (end-last).max(0)
        if dt.upper()>0:
            segments.append((dt,b,B,steps if i==0 or dt.upper()>2 else max(48,steps//128)))
        last = end
        b = [x+y for x,y in zip(b,w)]
        B = [[x+y for x,y in zip(r,s)] for r,s in zip(B,extra)]
    dt = (L-last).max(0)
    if dt.upper()>0:
        segments.append((dt,b,B,max(48,steps//128)))
    return segments


def crossing_enclosure(model, initial, layers, full, denominator, seed_L, seed_alpha, *, additional_mass=None, extra_thresholds=(), recenter=False):
    unit = arb("2.435e18")/denominator**2
    start = (unit/model["mz"]).log()
    table = [(start+l,w,B,f,s) for l,w,B,f,s in two["thresholds"](model,layers,full)]
    if additional_mass is not None:
        # 另一复弱双重态只计规范系数，额外 Yukawa 取零。
        indices = [arb(3)/10,arb(1)/2,arb(0)]
        casimirs = [arb(3)/20,arb(3)/4,arb(0)]
        extra = [[v*(4*casimirs[j]+(arb(2)*(0,2,3)[i]/3 if i==j else 0))
                  for j in range(3)] for i,v in enumerate(indices)]
        table.append(((additional_mass/model["mz"]).log(),model["higgs"],extra,0,1))
    for mass, one, matrix in extra_thresholds:
        table.append(((mass/model["mz"]).log(), one, matrix, 0, 1))
    if additional_mass is not None or extra_thresholds:
        table.sort(key=lambda row:float(row[0]))
        assert all(a[0]<b[0] for a,b in zip(table,table[1:]))
    L = arb(seed_L)+arb(0,".014")
    inverse = 1/arb(seed_alpha)+arb(0,".04")
    pi = arb.pi()
    assert all(0<l for l,*_ in table)
    even = [*layers[0],*(h for q in (2,4,6) for h in layers[q][3:])]

    def image(L,inverse,steps):
        # 以确定开启的阈值建立仿射主项，可能在盒内穿越的阈值保留正部。
        definite = [r for r in table if r[0]<L]
        possible = [r for r in table if not (r[0]>L)]
        segments = gauge_segments(L,possible,steps)
        corr,_ = two["flow"]([*initial[:2],inverse],segments,sm["D"],sm["C"],
                             min_inverse=7 if additional_mass is not None or recenter else 8)
        rhs, mat = [],[]
        for i in range(3):
            b0 = sm["BETA"][i]+sum(r[1][i] for r in definite)
            a = initial[i] if i!=2 else 0
            val = 2*pi*(a+corr[i])+sum(r[1][i]*r[0] for r in definite)
            val -= sum(w[i]*(L-l).max(0) for l,w,_,_,_ in possible[len(definite):])
            rhs.append([val])
            mat.append([b0,2*pi,-2*pi*int(i==2)])
        return arb_mat(mat).solve(arb_mat(rhs)).entries()
    if recenter:
        # 中点迭代只定位下一次区间检验的初盒。
        for _ in range(3):
            guess = image(L.mid(),inverse.mid(),1024)
            L = guess[0].mid()+arb(0,".014")
            inverse = guess[2].mid()+arb(0,".04")
    assert all(32*h*unit>model["mz"]*L.exp() for h in even)
    result = image(L,inverse,8192)
    assert L.contains_interior(result[0]) and inverse.contains_interior(result[2]),(result,L,inverse)
    for steps in (8192,16384,32768):
        L,inverse = L.intersection(result[0]),inverse.intersection(result[2])
        result = image(L,inverse,steps)
    L,inverse = L.intersection(result[0]),inverse.intersection(result[2])
    alpha = 1/inverse
    mu = model["mz"]*L.exp()
    print("十九标量" if full else "十六标量",denominator,"交点 / GeV 与反推强耦合",
          *[(v.lower().str(10),v.upper().str(10)) for v in (mu,alpha)])
    return mu, result[1], alpha


def solve(model, initial, layers, full, denominator, seed_L, seed_alpha):
    return crossing_enclosure(model,initial,layers,full,denominator,seed_L,seed_alpha)[2]


def verify():
    model,_,_,initial,_ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 192
    read = chiral["source"]["spectral_supplier"](model)
    rows = read(acb.zeta_zero(11661).imag)
    layers = {q: chiral["ordered"]([v for v,z in rows if z==q]) for q in (0,*model["charges"])}
    assert all(layers[q][2]<9<layers[q][3] for q in model["charges"])
    geometry = chiral["source"]["geometry"]
    clocks = geometry["phases"](acb.zeta_zero(11661).imag)
    assert set(clocks)=={724,804}
    for d in clocks:
        assert sum(k==d for _,_,k,_ in geometry["edges"])==6
    for full,n,L,alpha in ((True,724,"28.987",".1225"),(True,804,"28.703",".11655"),
                           (False,724,"29.040",".11838"),(False,804,"28.763",".11301")):
        out = solve(model,initial,layers,full,n,L,alpha)
        if not full and n==724:
            assert arb(".1183")<out<arb(".1185")
        else:
            assert out<arb(".117") or out>arb(".119")


if __name__=="__main__":
    verify()
