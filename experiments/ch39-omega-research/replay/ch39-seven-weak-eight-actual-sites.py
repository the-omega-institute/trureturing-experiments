"""同一八站局部酉可在第七至十三壳同时消去原生第三枢轴弱差。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    graph = lib("three-ray-chiral-boundary")
    native = lib("native-rh-ray-pvm")
    frame = lib("cyclic-two-characters-nine")["characters"](graph["geo"]["sites"])[0][0]
    floor, first, second = levels
    source = fifth*columns*acb_mat([[(first-floor).sqrt(),0],[0,(second-floor).sqrt()]])
    return verify(root,source,levels,graph,native,frame)


def verify(root, source, levels, graph, native, frame):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    gram = native["gram_cholesky"](root)
    rows, weights = [], []
    for number in range(2_306_888_100_092,2_306_888_100_099):
        _, basis = lib("adjacent-pvm-transport")["spectrum"](graph,6,number)
        rows.append(adj(basis))
        weights.append([v.real for v in native["readout"](basis,range(19),[14,8,7,6],gram[3],3)])
    ctx.prec = 256
    selected = acb_mat([m.tolist()[0] for m in rows])*frame
    inverse = selected.inv()
    centers = list(map(arb,(*lib("six-cyclic-weak-horizon")["Q"],"506.579672393")))
    width, planck, mz = arb(".00005"), arb("1.22089e19"), arb("91.1876")
    radii = [((q+arb(0,width.upper()))/planck-levels[0]).sqrt() for q in centers]
    fixed = (adj(frame)*source).tolist()[6]
    length = sum(abs(z)**2 for z in fixed).sqrt()
    direction = acb_mat([fixed])/length
    sides = [abs(inverse[6,j])*r for j,r in enumerate(radii)]
    left,right = sum(sides[:3]),sum(sides[3:])
    assert abs(left-right)<length<left+right
    x = (length**2+left**2-right**2)/(2*length)
    y = (left**2-x**2).sqrt()
    phases = [acb(x,y)/left if j<3 else acb(length-x,-y)/right for j in range(7)]
    target = acb_mat([[r*p*abs(inverse[6,j])/inverse[6,j]] for j,(r,p) in enumerate(zip(radii,phases,strict=True))])
    averages = inverse*target*direction
    support = ((3,4),(1,5),(8,),(9,),(10,),(11,))
    sites = [i for group in support for i in group]
    local = acb_mat([source.tolist()[i] for i in sites])
    prescribed = []
    for k,group in enumerate(support):
        active = [i for i in range(19) if abs(frame[i,k])>arb(".5")]
        assert set(group)<=set(active)
        scale = arb(len(active)).sqrt()
        prescribed.append([(scale*averages[k,h]-sum(source[i,h] for i in active if i not in group))
                           /arb(len(group)).sqrt() for h in (0,1)])
    average = acb_mat(prescribed)
    residual = adj(local)*local-adj(average)*average
    a = residual[0,0].real.sqrt()
    assert a>0
    b = residual[0,1]/a
    d2 = residual[1,1].real-abs(b)**2
    assert d2>0
    completion = acb_mat([[a,b],[0,d2.sqrt()]])
    # 两组双站各有一条 (1,-1)/√2 方向。补足两列的 Gram 后，局部等距可以延拓成八维酉。
    image = acb_mat(source)
    for k,group in enumerate(support):
        scale = arb(len(group)).sqrt()
        for t,i in enumerate(group):
            for h in (0,1):
                image[i,h] = (average[k,h]+((-1)**t)*completion[k,h])/scale if k<2 else average[k,h]

    reference = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    p0,p1 = (mz/planck)**2,(arb("1.176397e14")/planck)**2
    weak = native["weak_series"](56)
    large = lib("full-neutral-free-cutoff")["weak_column"]
    minimum, margin = arb(1), arb(1)
    assert len(rows)==len(weights)==len(centers)==7
    for row,w,q in zip(rows,weights,centers,strict=True):
        changed = row*image
        high_rows = [(levels[0]+sum(max(arb(0),abs(z).lower())**2 for z in r)).lower() for r in changed.tolist()[1:]]
        assert len(high_rows)==18 and all(z.is_finite() for z in high_rows)
        high = min(high_rows)
        assert high**2>p0
        error = (large(high**2,p0,p1)-reference).abs_upper()*sum(abs(v).upper() for v in w[1:])
        faces = [reference+w[0]*(weak(((q+sign*width)/mz)**2)-reference)+arb(0,error) for sign in (-1,1)]
        assert all(v.is_finite() for v in faces) and faces[0]>0>faces[1],(q,high,faces)
        minimum = min(minimum,high)
        margin = min(margin,faces[0].lower(),-faces[1].upper())
    print("八站支撑、Cholesky 两个主元平方下界、其余截止下界、七对弱差面的最小余量",
          sites,(a**2).lower().str(12),d2.lower().str(12),minimum.str(12),margin.str(12),flush=True)


if __name__=="__main__":
    root=Path(__file__).resolve().parent
    source=run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*source["verify"](root,return_operator="spectral"))
