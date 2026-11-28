"""有界标量移位下前 177 枚八站相消编号，以及跨三枚固定电流的残差。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx

SITES = (0, 2, 3, 5, 6, 9, 11, 13)
HITS = (120, 163, 177)
NORMALS = {36: (48, 88), 68: (71, -71), 83: (67, -75), 148: (59, -81)}


def remainder(green, radius):
    norm = sum(abs(z)**2 for z in green.entries()).sqrt().upper()
    column_sq = max(sum(abs(green[i,j])**2 for i in range(green.nrows())).upper() for j in SITES)
    assert radius*norm < 1
    return radius**2*column_sq*norm/(1-radius*norm)


def real_minima(green, derivative, radius, error):
    values = {True: [], False: []}
    for i in SITES:
        for j in SITES:
            special = {i,j} in ({6,3}, {6,5})
            lower = green[i,j].real-radius*abs(derivative[i,j].real)-error
            values[special].append(lower.lower())
    return min(values[False]), min(values[True])


def concentration(minima, threshold):
    good, bad = minima
    assert good > threshold and bad < good
    total = (good-threshold)/(good-bad)
    assert total > arb(1)/2
    return (total+(2*total-1).sqrt())/2


def enumerate_indices(laplacian):
    radius = arb(".02")
    interval = arb(0,radius.upper())
    hits, separated, direct = {}, [], []
    for index in range(1,178):
        green = laplacian(5,acb.zeta_zero(index).imag).inv()
        derivative = -green*green
        error = remainder(green,radius)
        minima = real_minima(green,derivative,radius,error)
        if min(minima) > 0:
            direct.append(min(minima))
            continue
        if index in NORMALS:
            assert concentration(minima, arb(0)) > arb(".8")
            x,y = NORMALS[index]
            normal_length = arb(x*x+y*y).sqrt()
            for i in SITES:
                for j in SITES:
                    z = (4*green[6,j]+green[i,j])/5
                    d = (4*derivative[6,j]+derivative[i,j])/5
                    assert x*z.real+y*z.imag > radius*abs(x*d.real+y*d.imag)+normal_length*error
            separated.append(index)
            continue
        assert index in HITS
        # 余项复模至多 error；矩形外包足以保持三个严格面积符号。
        noise = acb(arb(0,error.upper()),arb(0,error.upper()))
        x,y,z = [green[6,j]+interval*derivative[6,j]+noise for j in (0,3,5)]
        areas = [(u.conjugate()*v).imag for u,v in ((y,z),(z,x),(x,y))]
        assert all(t>0 for t in areas) or all(t<0 for t in areas)
        hits[index] = green, derivative
    assert len(direct)==170 and min(direct)>arb(".0000509")
    assert separated==list(NORMALS) and tuple(hits)==HITS
    print("全移位区间中的直接排除数、半平面排除与三个命中",len(direct),separated,tuple(hits))
    return hits


def vertices(matrix, concentration_bound, split):
    k = concentration_bound
    for i in SITES:
        row = {j:k*matrix[6,j]+(1-k)*matrix[i,j] for j in SITES}
        if split is None:
            yield from row.values()
            continue
        inner = split*row[3]+(1-split)*row[5]
        yield inner
        for j in SITES:
            if j not in (3,5):
                yield k*inner+(1-k)*row[j]


def fixed_pair(hits, first, second, radius, threshold, k, tests):
    green1, slope1 = hits[first]
    green2, slope2 = hits[second]
    errors = [remainder(g,radius) for g in (green1,green2)]
    # 假定两读数均不超过 threshold。在其中一个高度导出 a_6 与 b_3+b_5 的共同下界。
    close_bounds = [concentration(real_minima(g,d,radius,e),threshold)
                    for (g,d),e in zip((hits[first],hits[second]),errors,strict=True)]
    assert max(close_bounds)>k
    for splits, normal in tests:
        a,b,c,d = normal
        norm1,norm2 = arb(a*a+b*b).sqrt(),arb(c*c+d*d).sqrt()
        for split in splits:
            streams = [vertices(m,k,split) for m in (green1,green2,slope1,slope2)]
            for z,w,dz,dw in zip(*streams,strict=True):
                projection = a*z.real+b*z.imag+c*w.real+d*w.imag
                linear_loss = radius*(abs(a*dz.real+b*dz.imag)+abs(c*dw.real+d*dw.imag))
                assert projection-linear_loss-norm1*errors[0]-norm2*errors[1] > threshold*(norm1+norm2)
    print("编号对、独立移位半径和两个模长最大值的严格下界",first,second,radius.str(3),threshold.str(5))


def verify(root):
    laplacian = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))["laplacian"]
    ctx.prec = 192
    hits = enumerate_indices(laplacian)
    split = arb(".54912")
    tests = [((arb(0),split),(-336310,359998,368088,-349165)),
             ((split,arb(1)),(-336684,353531,368100,-355585))]
    fixed_pair(hits,120,163,arb(".02"),arb(".000005"),arb(".909"),tests)
    for idx,normal in ((120,(845,9060,-355,-10000)),(163,(1717,8931,-831,-10000))):
        for radius,threshold in ((".02",".00024"),(".15",".00004")):
            fixed_pair(hits,idx,177,arb(radius),arb(threshold),arb(".8"),[((None,),normal)])

    # 一般厄米扰动由预解式恒等式控制，静态下界也由同一顶点检验直接复核。
    for green,_ in hits.values():
        assert sum(abs(z)**2 for z in green.entries()).sqrt() < arb(".464")
        assert all(sum(abs(green[i,j])**2 for i in range(green.nrows())) < arb(".01583")
                   for j in SITES)
    fixed_pair(hits,120,163,arb(0),arb(".00003295"),arb(".909"),tests)
    for idx,normal in ((120,(845,9060,-355,-10000)),(163,(1717,8931,-831,-10000))):
        fixed_pair(hits,idx,177,arb(0),arb(".00024"),arb(".8"),[((None,),normal)])
    for size,gap in ((".002",".00003295"),(".015",".00024")):
        eta = arb(size)
        assert arb(".01583")*eta/(1-arb(".464")*eta) < arb(gap)


if __name__=="__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
