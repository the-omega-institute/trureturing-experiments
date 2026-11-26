"""两种正边权的 Higgs 动能与三低模夸克的顶底比。"""
import json
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adjoint(matrix):
    return matrix.conjugate().transpose()


def lowest(matrix, basis):
    block = adjoint(basis)*matrix*basis
    values, vectors = block.eig(right=True)
    index = min(range(len(values)),key=lambda j:float(values[j].real))
    level = values[index].real
    others = [v.real for j,v in enumerate(values) if j != index]
    assert all(level < x for x in others)
    vector = basis*acb_mat([[vectors[i,index]] for i in range(len(values))])
    return level, vector/sum(abs(a[0])**2 for a in vector.tolist()).sqrt(), others


def sorted_readings(left, right, profile, conjugate):
    h = profile.conjugate() if conjugate else profile
    entries = [abs(sum(a[i,0].conjugate()*h[i,0]*b[i,0] for i in range(19)))
               for a,b in zip(left,right)]
    entries.sort(key=float,reverse=True)
    assert all(x>y for x,y in zip(entries,entries[1:]))
    return entries


def spanning_tree(rotation):
    tree = set()
    for i,j in ((0,2),(1,14),(6,13),(9,2),(10,15),(11,2)):
        for _ in range(3):
            tree.add(tuple(sorted((i,j))))
            i,j = rotation[i],rotation[j]
    return tree


def verify():
    ctx.prec = 224
    root = Path(__file__).parent
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    common = run_path(str(root/"ch39-all-common-shell-ranks.py"))
    ctx.prec = 224
    sites = shell["geo"]["sites"]
    rotation = [sites.index((y,z,x)) for x,y,z in sites]
    bases = list(common["character_bases"](shell))
    height = acb.zeta_zero(11431).imag
    matter = {}
    for q in (1,4,-2):
        matter[q] = []
        lap = shell["laplacian"](q,height)
        for basis in bases:
            level, mode, rest = lowest(lap,basis)
            assert level<9 and all(x>9 for x in rest)
            matter[q].append(mode)

    edges = {(i,j) for i,j,_,_ in shell["geo"]["edges"]}
    tree = spanning_tree(rotation)
    assert len(tree) == 18 and len(edges) == 135 and tree <= edges
    connected = {14}
    for _ in range(18):
        connected.update(j for i,j in tree if i in connected)
        connected.update(i for i,j in tree if j in connected)
    assert len(connected) == 19
    original = shell["laplacian"](3,height)
    base = acb_mat(19,19)
    degree = [arb(0) for _ in sites]
    for i,j in edges:
        base[i,j] = original[i,j] / (1 if (i,j) in tree else 1000)
        base[j,i] = base[i,j].conjugate()
        degree[i] += abs(base[i,j])
        degree[j] += abs(base[i,j])
    # 边模长与荷无关，故 Gershgorin 界适用于共用重权的五种荷。
    assert all(2*x < arb("8.013") for x in degree)
    shared_modes = {}
    for charge in (1,4,-2,3,6):
        raw = shell["laplacian"](charge,height)
        weighted = acb_mat(19,19)
        for i in range(19):
            weighted[i,i] = degree[i]
        for i,j in edges:
            weighted[i,j] = raw[i,j] / (1 if (i,j) in tree else 1000)
            weighted[j,i] = weighted[i,j].conjugate()
        shared_modes[charge] = []
        for basis in bases:
            level,mode,rest = lowest(weighted,basis)
            assert level < arb(".123") and all(v>arb(".322") for v in rest)
            shared_modes[charge].append(mode)
        if charge == 3:
            kinetic = weighted

    projectors = {q: sum((v*adjoint(v) for v in shared_modes[q]),acb_mat(19,19))
                  for q in (1,4,-2)}
    case = json.loads((root/"ch39-common-tree-certificate.json").read_text())
    bound = case["bound"]
    assert bound == 3
    gram = run_path(str(root/"ch39-higgs-phase-vacuum.py"))["gram"](
        rotation,{"gram":[arb(n)/10 for n in case["gram"]]})
    positive = run_path(str(root/"ch39-three-mode-higgs-hierarchy.py"))["positive"]
    positive(gram)
    v = [arb(n)/10**7 for n in case["first"]]
    certificate = acb_mat(19,19)
    for i in range(19):
        for j in range(19):
            factor = gram[i,i]+gram[j,j]-2*gram[i,j]+acb(0,v[i]-v[j])
            kd = projectors[1][i,j]*projectors[-2][i,j].conjugate()
            ku = projectors[1][i,j].conjugate()*projectors[4][i,j]
            certificate[i,j] = factor*kinetic[i,j]+300*kd-900*ku/bound**2
    positive(certificate)
    for i,value in zip((0,1,6,9,10,11,14),(1392,29310,540,1743,2472,458,144)):
        for _ in range(3):
            base[i,i] = arb(value)/1000
            i = rotation[i]

    reference = arb(8573)/141
    def read(t):
        matrix = acb_mat(base)
        matrix[14,14] -= t
        level,h,others = lowest(matrix,bases[0])
        for basis in bases[1:]:
            others.extend(v.real for v in (adjoint(basis)*matrix*basis).eig(multiple=True))
        assert level < 0 and all(x-level>arb(".007") for x in others)
        up = sorted_readings(matter[1],matter[4],h,True)
        down = sorted_readings(matter[1],matter[-2],h,False)
        return up,down

    lo,hi = arb(38012548)/10**11,arb(38012550)/10**11
    left,ld = read(lo)
    right,rd = read(hi)
    assert ld[0]>0 and rd[0]>0 and left[0]/ld[0]<reference<right[0]/rd[0]
    # Weyl 和低线与正交补的投影界，允许整个闭区间与左端比较。
    e = hi-lo
    assert arb(".007")-2*e>0
    mode_change = 2*e/(arb(".007")-e)
    assert ld[0]-mode_change > arb(".00018")
    assert (left[1]-mode_change)/(left[0]+mode_change) > arb(".536")
    assert (left[1]+mode_change)/(left[0]-mode_change) < arb(".538")
    up,down = read(arb(515)/10**6)
    assert down[0]>0 and 4700<up[0]/down[0]<4900
    print("共用边权的五荷共同截断与全部非零局部真空顶底比上界", "1/5", bound)
    print("保留与减重边数、参考顶底比的极点调节区间、固定有理调节的顶底比",
          len(tree),len(edges)-len(tree),lo.str(12),hi.str(12),(up[0]/down[0]).str(12))


if __name__ == "__main__":
    verify()
