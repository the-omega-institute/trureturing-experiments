"""从联合零点十六重层的原始根核验非钟边、非辐边和正向电导改动的必要下界。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


def support_exclusion(p, v, s, edges, vertices, outside, count):
    # 若专用边外至多 count 条超过 1，对正向增量和 u 的严格凸下界。
    h = arb(19)/(2*vertices*(19-vertices))
    t = (19-2*vertices)*s/19-outside
    minimum_at = (p/edges+h*t)/(1/arb(edges)+1/arb(count)+h)
    assert 0 < minimum_at < p and minimum_at < t
    return ((p-minimum_at)**2/edges+minimum_at**2/count+h*(t-minimum_at)**2-v).lower()


def verify(root):
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 256
    h = acb.zeta_zero(2306888100085).imag
    delta = arb(23)/25*arb(2)**-20/(h*h+shell["soft"]["point"]["study"]["offset"])
    c = 15-delta/3
    edges = [(i,j) for i,j,_,_ in shell["geo"]["edges"]]
    clock = [(i,j) for i,j,d,_ in shell["geo"]["edges"] if d in (724,804)]
    pole = shell["soft"]["pole"]
    spokes = [e for e in edges if pole in e]
    three = {6,7,8}
    four = three | {pole}
    assert len(edges)==135 and len(clock)==12 and len(spokes)==18 and not set(clock)&set(spokes)
    assert {i for e in clock for i in e}==set(range(9))
    assert all(set(e)&three for e in clock) and all(set(e)&four for e in clock+spokes)
    assert [sum(not set(e)&core for e in edges) for core in (three,four)] == [84,69]
    neighbors = [[j for e in edges if i in e for j in e if j!=i] for i in range(19)]
    deg = [len(n) for n in neighbors]
    def signless(x):
        return [deg[i]*x[i]+sum(x[j] for j in neighbors[i]) for i in range(19)]
    perron = [1]*19
    for _ in range(32):
        perron = signless(perron)
    b = max(arb(t)/x for t,x in zip(signless(perron),perron))
    assert b<arb("30.386")
    tr = sum(deg)
    squares = sum(n*n for n in deg)+tr
    lower_three = arb(3*tr)/19-(arb(3*16)/19*(squares-arb(tr*tr)/19)).sqrt()
    assert (tr,squares)==(270,4290) and lower_three>arb("8.79639")
    max_exterior_five = max(sum(not (i in C or j in C) for i,j in edges)
                            for C in combinations(range(19),5))
    for q, counts_fixed, counts_free in ((2,(25,15,8),(24,14,7)),
                                        (4,(25,14,8),(24,14,7)),
                                        (6,(25,14,7),(24,14,7))):
        lap = shell["laplacian"](q,h)
        values = lap.eig(multiple=True)
        assert all(z.real>0 and z.imag.contains(0) for z in values)
        raw = sorted((z.real for z in values),key=float)
        masses = sorted((x*abs(x-c) for x in raw),key=float)
        assert all(x<y for x,y in zip(masses,masses[1:]))
        assert 32*masses[3]>c*c/4
        heavy = [(c+(c*c+128*m).sqrt())/2 for m in masses[3:]]
        s = sum(heavy)/2
        assert 2*s/19>heavy[0] and heavy[0]>b
        variance = sum(y*y for y in heavy)-(2*s)**2/19
        v = variance/2-2*s+135
        # 将三条原软边补至 1 后，增加矩阵范数为 4δ。
        base_top = raw[-1]+4*delta
        base_16 = sum(raw[-16:])+16*4*delta
        outside_gain = s-sum(heavy[-4:])-69
        assert outside_gain > 80 and s>sum(heavy[-5:])+max_exterior_five
        thresholds=(heavy[0]/base_top,heavy[0]/b,(s-sum(heavy[-4:]))/69)
        print("荷、原谱顶、最小重根、最大重根、四点外净增量下界、原相/换相16条边门槛及四点外单边门槛",
              q,*[val.str(13) for val in (raw[-1],heavy[0],heavy[-1],outside_gain,*thresholds)])
        for name,upper,counts in (("原相位",base_16,counts_fixed),("任意换相",tr-lower_three,counts_free)):
            p = s-upper/2
            minors = (p*p/(counts[0]-1)-v,
                      support_exclusion(p,v,s,12,3,84,counts[1]-1),
                      support_exclusion(p,v,s,30,4,69,counts[2]-1))
            assert all(x>0 for x in minors)
            print("相位、超过单位电导的总边数/非钟边数/非钟非辐边数、正向总量与三项严格余量",
                  name,counts,p.str(13),*[x.str(12) for x in minors])
    print("单位权任意相位的谱顶上界、底三根和下界、五点外原边数最大值",
          b.str(15),lower_three.str(14),max_exterior_five)


if __name__=="__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
