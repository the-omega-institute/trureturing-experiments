"""标准零点谱权乘以一条壳辐边上的一般正读出。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


root = Path(__file__).parent
scale = run_path(str(root/"ch39-scale-and-zero-weights.py"))
shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))


def compression(charge, height, sites):
    lap = shell["laplacian"](charge,height)
    center = sum(lap[i,i].real for i in range(19))/18
    values,vectors = lap.eig(right=True,algorithm="rump")
    assert all(not z.overlaps(w) for j,z in enumerate(values) for w in values[:j])
    pole, tip = sites
    a,b,c = arb(0),arb(0),acb(0)
    for j,z in enumerate(values):
        mass = abs(z.real*(z.real-center))
        assert mass<9 or mass>9
        if mass<9:
            norm = sum(abs(vectors[i,j])**2 for i in range(19))
            u,v = vectors[pole,j],vectors[tip,j]
            a += abs(u)**2/norm
            b += abs(v)**2/norm
            c += u*v.conjugate()/norm
    return [a,b,2*c.real,2*c.imag]


def sector_bound(q, witness, pair, heights, ratios, total, sites):
    support = sorted({1,2,3,4,5,witness})
    rows = {n:compression(q,heights[n-1],sites) for n in support}
    a,b,c,d = rows[witness]
    eigen_floor = (a*b-(c*c+d*d)/4)/(a+b)
    assert 0 < eigen_floor < 1
    base_weight = 1/(total*heights[witness-1]**2)
    assert base_weight*eigen_floor > arb("1e-9")
    trace_cap = 1/(base_weight*eigen_floor)
    if pair:
        i,j = pair
        pole_error = {k:(abs(rows[k][0])+abs(rows[k][2])+abs(rows[k][3]))*trace_cap for k in pair}
        diff = rows[j][1]-rows[i][1]
        assert not diff.contains(0)
        tip_cap = (abs(ratios[j]-ratios[i])+pole_error[i]+pole_error[j])/abs(diff)
        fraction = 1-ratios[i]+pole_error[i]+rows[i][1]*tip_cap
    else:
        differences = arb_mat([[b-a for a,b in zip(rows[1],rows[n],strict=True)] for n in range(2,6)])
        response = differences.solve(arb_mat([[ratios[n]-ratios[1]] for n in range(2,6)]))
        fraction = 1-ratios[1]+sum(value*response[j,0] for j,value in enumerate(rows[1]))
    return fraction.upper()


def verify():
    ctx.prec,ctx.cap=1536,3
    total=scale["local"]["mask"]["time_data"]()[0]
    times=scale["clock"]["times"]()
    windows=scale["local"]["required_windows"](times)
    weights=scale["small_scale_weights"](times,windows[0][1],arb(2)**-176)
    heights=[z.imag for z in acb.zeta_zeros(1,35)]
    ratios={i:w*total*heights[i-1]**2 for i,w in enumerate(weights,1)}
    ctx.prec=320
    points=shell["geo"]["sites"]
    sites=tuple(points.index(v) for v in ((17,17,17),(13,23,13)))
    cases = ((1,1,None,"6.01e-9"), (2,1,(2,3),"2.73e-14"), (3,1,None,"4.87e-8"),
             (4,4,(1,2),"3.24e-14"), (6,35,(1,2),"4.16e-13"))
    for q,witness,pair,ceiling in cases:
        upper=sector_bound(q,witness,pair,heights,ratios,total,sites)
        print("荷、固定辐边任意正探针的总权份额上界",q,upper.str(13))
        assert upper<arb(ceiling)


if __name__=="__main__":
    verify()
