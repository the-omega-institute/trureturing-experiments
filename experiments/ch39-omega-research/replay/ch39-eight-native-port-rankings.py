"""用四个固定端口及原生射线 Gram 权核验八站截取。"""
from collections import Counter
from itertools import combinations
from math import comb
from pathlib import Path
from runpy import run_path

from flint import arb, ctx, fmpq, fmpq_mat


def tied_prefix(blocks, winners):
    """整列并列同时保留，报告全部合法的并列截取与其中的达标支撑数。"""
    strict = []
    for group in blocks:
        take = 8-len(strict)
        if len(group) >= take:
            selections = {tuple(sorted(strict+list(xs))) for xs in combinations(group,take)}
            assert len(selections) == comb(len(group),take)
            return len(selections), len(selections & winners)
        strict.extend(group)
    raise ValueError("站数不足")


def port_fields(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    shell = lib("three-ray-chiral-boundary")
    sites = shell["geo"]["sites"]
    ports = (shell["soft"]["pole"], *shell["boundary_sites"]())
    assert ports == (14,8,7,6) and len(sites) == 19
    lap = fmpq_mat(19,19)
    for i,j,_,_ in shell["geo"]["edges"]:
        lap[i,i] += 1
        lap[j,j] += 1
        lap[i,j] -= 1
        lap[j,i] -= 1
    mean = fmpq_mat([[fmpq(1,19)]*19 for _ in sites])
    green = (lap+mean).inv()
    second = lap*lap

    def dipole(matrix,i,j):
        return (matrix[i,i]+matrix[j,j]-matrix[i,j]-matrix[j,i])/2

    features = {}
    features["平方弦长"] = [tuple(fmpq(sum((x-y)**2 for x,y in zip(sites[i],sites[j],strict=True))) for j in ports)
                             for i in range(19)]
    for name,m in (("一阶热矩",lap),("二阶热矩",second),("半有效电阻",green)):
        features[name] = [tuple(dipole(m,i,j) for j in ports) for i in range(19)]
    ctx.prec = 256
    gram = lib("native-rh-ray-pvm")["gram_cholesky"](root)
    born = [w**2 for w in gram[3]]
    assert all(w>0 for w in born) and sum(born).contains(1)
    return sites,ports,features,born


def interior_witness(support, sites):
    """在所选线段上找一个严格位于补集四面体内部的有理点。"""
    outside = [i for i in range(19) if i not in support]
    tetrahedra = []
    for tetra in combinations(outside,4):
        matrix = fmpq_mat([(*sites[i],1) for i in tetra]).transpose()
        if matrix.det():
            tetrahedra.append((tetra,matrix.inv()))
    pairs = list(combinations(support,2))
    for numerator,denominator in ((1,2),(1,3),(2,3)):
        columns = [[fmpq((denominator-numerator)*sites[i][h]+numerator*sites[j][h],denominator)
                    for h in range(3)]+[1] for i,j in pairs]
        points = fmpq_mat(columns).transpose()
        for tetra,inverse in tetrahedra:
            coords = inverse*points
            for k,pair in enumerate(pairs):
                weights = tuple(coords[h,k] for h in range(4))
                if min(weights)>0:
                    assert sum(weights)==1
                    assert all(sum(w*sites[i][h] for w,i in zip(weights,tetra,strict=True))==columns[k][h]
                               for h in range(3))
                    return pair,numerator,denominator,tetra,weights
    raise AssertionError(("无证书",support))


def geometric_obstructions(winners, sites, features):
    assert all(sum(a*a for a in x)==867 for x in sites)
    pole = (17,17,17)
    dot = lambda p,q: sum(a*b for a,b in zip(p,q,strict=True))
    # ch39-arc-relative-determinant.py 的原始球面三角相位，分母在全壳上严格为正。
    denominators = [867+dot(pole,x)+dot(x,y)+dot(y,pole) for x in sites for y in sites]
    assert min(denominators)>0
    print("原始球面三角相位分母括号下界",min(denominators),flush=True)
    counts=Counter()
    example=(1,3,4,5,8,9,10,11)
    for support in sorted(winners):
        pair,num,den,tetra,bary=interior_witness(support,sites)
        counts[den]+=1
        if support==example:
            print("示例线段、分点、补集四面体、正重心坐标",pair,num,den,tetra,bary,flush=True)
        outside=[i for i in range(19) if i not in support]
        for name,rows in features.items():
            if name=="平方弦长":
                continue
            for sign in (-1,1):
                # 一个内外站对在全部四端口都严格逆序，任意非零非负混合因此失败。
                assert any(all(sign*(a-b)>0 for a,b in zip(rows[i],rows[j],strict=True))
                           for i in support for j in outside),(support,name,sign)
    assert counts=={2:29,3:3}
    print("非恒定仿射排序的中点、三分点排除数",counts[2],counts[3],flush=True)
    print("每种无磁图端口量、每一排序方向均有四端口严格支配反例",len(winners),flush=True)


def single_spectral_lines(root, winners):
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    sites = shell["geo"]["sites"]
    cyclic = [sites.index((b,c,a)) for a,b,c in sites]
    orbits = sorted({tuple(sorted({i,cyclic[i],cyclic[cyclic[i]]})) for i in range(19)})
    partial = [group for group in orbits if all(0<len(set(group)&set(s))<len(group) for s in winners)]
    assert len(partial)==5
    number = run_path(str(root/"ch39-fifth-flavor-cutoff-chart.py"))["NUMBER"]
    spectrum = run_path(str(root/"ch39-adjacent-pvm-transport.py"))["spectrum"]
    _,basis = spectrum(shell,6,number)
    separation = arb(1)
    for line in range(19):
        values=[abs(basis[group[0],line])**2 for group in partial]
        assert all(v.is_finite() for v in values)
        gap=max(v.lower() for v in values)-min(v.upper() for v in values)
        assert gap>arb(".04266")
        separation=min(separation,gap)
    print("第五壳十九条单谱线的五轨道代表 Born 权，最小极差",separation.str(12),flush=True)


def audit(winners,sites,ports,features,born):
    geometric_obstructions(winners,sites,features)
    evidence = {}
    for name,features_by_site in features.items():
        for anchor in range(5):
            signatures = {}
            for site,row in enumerate(features_by_site):
                signature = (row[anchor],) if anchor < 4 else row
                signatures.setdefault(signature,[]).append(site)
            def score(key):
                if anchor < 4:
                    return arb(key[0])
                return sum(w*arb(x) for w,x in zip(born,key,strict=True))
            keys = sorted(signatures,key=lambda k:score(k).mid())
            assert all(score(a)<score(b) for a,b in zip(keys,keys[1:]))
            blocks = [signatures[k] for k in keys]
            label = str(ports[anchor]) if anchor<4 else "第三源Born期望"
            low,high = tied_prefix(blocks,winners),tied_prefix(blocks[::-1],winners)
            assert low[1] == high[1] == 0,(name,label,low,high)
            evidence[name,label] = low,high
            print(name,label,"低八、高八各自的候选数和可行数",low,high,flush=True)
    # 四个具名端口保持不动，十五个未带标签的站位任意置换。
    port_set = set(ports)
    hits = Counter(tuple(sorted(set(w)&port_set)) for w in winners)
    assert hits == {(6,):6,(7,):9,(8,):17} and comb(15,7)==6435
    rotation = [sites.index((b,c,a)) for a,b,c in sites]
    def rotate(s):
        return tuple(sorted(rotation[i] for i in s))
    representatives = {min(w,rotate(w),rotate(rotate(w))) for w in winners}
    occupancy = Counter(sum(s in winners for s in (w,rotate(w),rotate(rotate(w)))) for w in representatives)
    assert occupancy=={1:8,2:12}
    print("三循环不变的八站抽样命中上界 2/3；仅标记四个 Gram 端口的抽样命中上界 17/6435",flush=True)
    return evidence


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    fields = port_fields(root)
    source = run_path(str(root/"ch39-fifth-full-positive-joint.py"))["verify"](root,return_operator="spectral")
    winners = run_path(str(root/"ch39-eight-optimal-classification.py"))["experiment"](root,*source)
    audit(winners,*fields)
    single_spectral_lines(root,winners)
