"""固定同一编号零点、原边身份与两个整数标尺，在有理球面变形上重算五荷谱、两圈交点及四层曲率。"""
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx, fmpz_mat


NUMBER = 2_306_888_100_085
POLE = (17, 17, 17)


def enclose(v):
    return arb(v.numerator)/v.denominator


def links(graph, sites, phases, softened, deficit):
    pole = graph["sites"].index(POLE)
    dot, det = graph["dot"], graph["det"]
    prepared = []
    for i, j, original_length, _ in graph["edges"]:
        x, y = sites[i], sites[j]
        re = 867*(867+dot(POLE,x)+dot(x,y)+dot(y,POLE))
        im = -det(POLE,x,y)
        z = acb(enclose(re), arb(867).sqrt()*enclose(im))
        signed = (-1 if original_length == 200 else 1)*z/abs(z)
        unit = signed*phases.get(original_length,1)
        if pole in (i,j):
            assert (unit-1).contains(0)
            unit = acb(1)
        prepared.append((i,j,unit,1-deficit if (i,j) in softened else arb(1)))
    return prepared


def laplacian(prepared,q):
    matrix = acb_mat(19,19)
    for i,j,unit,weight in prepared:
        matrix[i,i] += weight
        matrix[j,j] += weight
        entry = -weight*unit**q
        matrix[i,j],matrix[j,i] = entry,entry.conjugate()
    return matrix


def neutral_levels(edges,deficit):
    flat = fmpz_mat(19,19)
    for i,j in edges:
        flat[i,i] += 1
        flat[j,j] += 1
        flat[i,j],flat[j,i] = -1,-1
    spectrum = flat.charpoly().complex_roots()
    assert sum(m for z,m in spectrum if z==0)==1
    # 去掉常数零模。三辐边星图的拉普拉斯范数为 4，Weyl 界给出所有余下质量。
    loss = arb(0).union(4*deficit)
    return [z.real-loss for z,m in spectrum if z!=0 for _ in range(m)]


def verify(root):
    physical = run_path(str(root/"ch39-first-three-mode-observables.py"))
    model,initial = physical["model_inputs"]()
    chiral = physical["joint"]["chiral"]
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    shape = run_path(str(root/"ch39-metric-pair-flexibility.py"))
    graph = shell["geo"]
    points = graph["sites"]
    read = run_path(str(root/"ch39-scale-free-geometric-pair.py"))["readout"]
    geometry = read(points)
    edges = geometry[1]
    assert geometry[0]==(32,724) and geometry[2]==(32,724,804)
    assert edges=={(i,j) for i,j,_,_ in graph["edges"]}
    ctx.prec,ctx.cap = 256,3
    height = acb.zeta_zero(NUMBER).imag
    phases = graph["phases"](height)
    assert set(phases)=={724,804}
    soft = shell["soft"]
    pole,tuned = soft["pole"],soft["star"]["orbit"]
    softened = {(i,j) for i,j in edges if pole in (i,j) and (i in tuned or j in tuned)}
    assert len(softened)==3
    deficit = arb(23)/25*arb(2)**-20/(height**2+soft["point"]["study"]["offset"])
    center = arb(15)-deficit/3
    cut = 9+(1-deficit)**2
    zero = neutral_levels(edges,deficit)
    assert len(zero)==18 and all(x>0 for x in zero)
    pairs = list(combinations(range(19),2))

    def exact_classes(sites):
        ds = {(i,j):shape["distance"](sites[i],sites[j]) for i,j in pairs}
        return tuple({e for e in pairs if ds[e]==d} for d in (244,724,804))

    classes = exact_classes(points)
    old_layer = None
    print("实际正零点序号、高度",NUMBER,height.str(24),flush=True)
    for denominator in (None,200,25,-25):
        parameter = Fraction(1,denominator) if denominator is not None else Fraction(0)
        sites = shape["deform"](points,parameter,parameter)
        assert read(sites)==geometry and exact_classes(sites)==classes
        assert sum(x!=y for x,y in zip(sites,points))==(6 if parameter else 0)
        assert all(graph["dot"](p,p)==867 for p in sites)
        prepared = links(graph,sites,phases,softened,deficit)
        levels = {0:zero}
        for q in model["charges"]:
            matrix = laplacian(prepared,q)
            if parameter==0:
                native = shell["laplacian"](q,height,phases=phases)
                assert all(v.contains(0) for v in (matrix-native).entries())
            eig = matrix.eig(multiple=True)
            assert len(eig)==19 and all(x.imag.contains(0) and x.real>0 for x in eig)
            levels[q] = [x.real for x in eig]
        layers = {q:chiral["ordered"]([x*abs(x-center) for x in row]) for q,row in levels.items()}
        print("有理参数 a=b、变形后六站的有序代表",parameter,sites[points.index((7,17,23))])
        for q in model["charges"]:
            third,fourth = layers[q][2:4]
            print("荷、第三与第四层", q, third.str(16), fourth.str(16))
        if denominator==-25:
            counts = []
            for q in model["charges"]:
                assert all(x<cut or cut<x for x in layers[q])
                counts.append(sum(x<cut for x in layers[q]))
            assert max(layers[q][2].lower() for q in model["charges"])>min(layers[q][3].upper() for q in model["charges"])
            assert list(model["charges"])==[1,4,2,3,6] and counts==[3,3,3,4,2]
            print("最小割下五荷的实际轻模个数",counts)
            continue
        assert all(layers[q][2]<9<cut<10<layers[q][3] for q in model["charges"])
        mass = layers[2][15]
        energy,alpha,curve = physical["layer_observables"](model,initial,layers)
        low,high = alpha.lower(),alpha.upper()
        print("荷二第十六层、两圈强耦合两端、四层曲率比",mass.str(18),low.str(14),high.str(14),curve.str(18),flush=True)
        print("两圈交点能标 / GeV",energy.lower().str(14),energy.upper().str(14))
        low_mass = (arb("246.22")**2/(arb("1.22089e19")*32*arb("2.435e18")/724**2*mass))**2
        coefficient = (9*arb(3).log()-12*arb(2).log())/2
        assert low_mass/(arb(".03")*coefficient)<arb("1e-57")
        assert arb(".117")<alpha<arb(".119")
        if parameter==0:
            old_layer = mass,curve,alpha
            assert arb(".11837517")<alpha<arb(".11838132")
            assert arb(".99999336")<curve<arb(".99999337")
        else:
            m0,c0,original_alpha = old_layer
            # 归一化低层 a 正比于质量的负二次幂，四层主项正比于 a²。
            assert (curve/c0).overlaps((m0/mass)**4)
            if denominator==200:
                assert arb(".11837265")<alpha<arb(".11837882")
                assert arb("1.000368")<curve<arb("1.000369")
            if denominator==25:
                assert arb(".11836369")<alpha<arb(".11837099") and alpha<original_alpha
                assert arb("3.7294e14")<energy<arb("3.7304e14")
                assert arb("1.002513")<curve<arb("1.002514")


if __name__=="__main__":
    root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)
