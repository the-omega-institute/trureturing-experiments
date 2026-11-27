"""平方弦长 32 匹配与 52 对称双叶星上的同一输运，同时匹配第七、第八壳。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat


def prepare(root, lib, fifth, columns, levels):
    sites = lib("three-ray-chiral-boundary")["geo"]["sites"]
    edges = {32: [], 52: []}
    for i,j in combinations(range(19),2):
        dist = sum((a-b)**2 for a,b in zip(sites[i],sites[j],strict=True))
        if dist in edges:
            edges[dist].append((i,j))
    assert edges[32] == [(0,2),(1,3),(4,5)]
    assert edges[52] == [(6,9),(6,13),(7,10),(7,15),(8,16),(8,17)]
    assert set(sum(edges[32],())).isdisjoint(sum(edges[52],()))
    geometries = []
    for number in (2_306_888_100_092, 2_306_888_100_093):
        base,_,_,basis,station,weights,heat,field = lib("fifth-seventh-shortest-site-gates")["geometry"](
            root, fifth, columns, levels, number=number)
        def rows(matrix, indices):
            return acb_mat([[matrix[i,j] for j in range(matrix.ncols())] for i in indices])
        def adj(m):
            return m.conjugate().transpose()
        def star_vectors(mat):
            middle = rows(mat, (6,7,8))
            bright = (rows(mat,(9,10,16)) + rows(mat,(13,15,17)))/arb(2).sqrt()
            return middle, bright
        bv,bw = star_vectors(basis)
        sv,sw = star_vectors(station)
        star = (adj(bv)*sv, adj(bw)*sw, adj(bw)*sv, -adj(bv)*sw)
        ends = tuple(zip(*edges[32],strict=True))
        match = []
        for sign in (1,-1):
            v = rows(basis,ends[0]) + sign*rows(basis,ends[1])
            s = rows(station,ends[0]) + sign*rows(station,ends[1])
            match.append(adj(v)*s/2)
        fixed = base - star[0]-star[1]-match[0]-match[1]
        geometries.append((fixed,star,match,weights,heat,field))
    return geometries


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    geometry = prepare(root,lib,fifth,columns,levels)
    gravity = lib("shortest-common-complex-gate")["gravity_system"]
    beta, gamma = arb("-2.375"), arb("-.125")
    imaginary = acb(0,1)
    def system(point, values):
        theta,alpha,angle,phase = point
        s,c = theta.sin_cos()
        ph = [(imaginary*t).exp() for t in (gamma+alpha,gamma-alpha,gamma+beta,gamma-beta)]
        v = [k*t for k,t in zip(ph,(c,c,s,s),strict=True)]
        dt = [k*t for k,t in zip(ph,(-s,-s,c,c),strict=True)]
        da = [imaginary*v[0],-imaginary*v[1]]
        plus,minus = [(imaginary*(phase+t)).exp() for t in (angle,-angle)]
        fs,js = [],[]
        for fixed,blocks,match,w,heat,field in geometry:
            p = fixed+sum((a*b for a,b in zip(v,blocks,strict=True)),0*fixed)
            p += plus*match[0]+minus*match[1]
            der = [sum((a*b for a,b in zip(dt,blocks,strict=True)),0*fixed),
                   da[0]*blocks[0]+da[1]*blocks[1], imaginary*(plus*match[0]-minus*match[1]),
                   imaginary*(plus*match[0]+minus*match[1])]
            f,j = gravity(p,der,levels,w,heat,field,values)
            fs.extend(f.tolist())
            js.extend(j.tolist())
        return arb_mat(fs),arb_mat(js)
    center = list(map(arb,("2.56245335564450",".25645082888507","1.90894913731359",".97134727724686")))
    values, jac = system(center,True)
    pre = jac.mid().inv()
    eye = arb_mat([[int(i==j) for j in range(4)] for i in range(4)])
    move = max(v.abs_upper() for v in (pre*values).entries())
    radius = arb(".0002")
    box = [v+arb(0,radius) for v in center]
    _,wide = system(box,False)
    defect = max(sum(abs(a).upper() for a in row) for row in (eye-pre*wide).tolist())
    ratio = move/radius+defect
    assert not pre.det().contains(0) and ratio < arb(".795")
    print("四式的共同半径、预条件中心偏移、导数误差、自映比",
          *[v.str(15) for v in (radius,move,defect,ratio)],flush=True)
    return box,ratio


def verify(root):
    source = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    return experiment(root,*source["verify"](root,return_operator="spectral"))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
