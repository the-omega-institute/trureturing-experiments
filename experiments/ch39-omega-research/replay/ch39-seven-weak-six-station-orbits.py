"""在三种六轨道支撑上，各有使第七至十三壳原生弱差同时为零的一次酉输运。"""
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
    frames, _ = lib("cyclic-two-characters-nine")["characters"](graph["geo"]["sites"])
    gram = native["gram_cholesky"](root)
    floor, first, second = levels
    weighted = fifth*columns*acb_mat([[(first-floor).sqrt(),0], [0,(second-floor).sqrt()]])
    rows, weights = [], []
    for number in range(2_306_888_100_092,2_306_888_100_099):
        _, basis = lib("adjacent-pvm-transport")["spectrum"](graph, 6, number)
        rows.append(adj(basis))
        weights.append([v.real for v in native["readout"](basis, range(19), [14,8,7,6], gram[3], 3)])
    ctx.prec=256
    pivot = frames[0]
    B = acb_mat([row.tolist()[0] for row in rows])*pivot
    inv = B.inv()
    centers = list(map(arb, (*lib("six-cyclic-weak-horizon")["Q"], "506.579672393")))
    width, planck, mz = arb(".00005"), arb("1.22089e19"), arb("91.1876")
    radii = [((q+arb(0,width.upper()))/planck-floor).sqrt() for q in centers]
    weak = native["weak_series"](56)
    large = lib("full-neutral-free-cutoff")["weak_column"]
    reference = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    p0,p1 = (mz/planck)**2,(arb("1.176397e14")/planck)**2
    source = adj(pivot)*weighted
    for absent, group in ((3, (2,3)), (4, (2,3)), (6, (0,1,2))):
        untouched = [i for i in range(19) if abs(pivot[i,absent])>arb(".5")]
        active = [i for i in range(19) if i not in untouched]
        s = sum(abs(z)**2 for z in source.tolist()[absent]).sqrt()
        direction = acb_mat([source.tolist()[absent]])/s
        lengths = [abs(inv[absent,j])*radii[j] for j in range(7)]
        left = sum(lengths[j] for j in group)
        right = sum(z for j,z in enumerate(lengths) if j not in group)
        assert abs(left-right)<s<left+right
        x = (s**2+left**2-right**2)/(2*s)
        y = (left**2-x**2).sqrt()
        phases = [acb(x,y)/left if j in group else acb(s-x,-y)/right for j in range(7)]
        target = acb_mat([[radii[j]*phases[j]*abs(inv[absent,j])/inv[absent,j]] for j in range(7)])
        trivial = inv*target*direction
        local_trivial = acb_mat([line for j,line in enumerate(trivial.tolist()) if j!=absent])
        local_source = acb_mat([weighted.tolist()[i] for i in active])
        remaining = adj(local_source)*local_source - adj(local_trivial)*local_trivial
        a = remaining[0,0].real.sqrt()
        assert a>0
        b = remaining[0,1]/a
        d2 = remaining[1,1].real-abs(b)**2
        assert d2>0
        # 第一条三站轨道与三种不动轨道均不相交，分别取该轨道的两种非平凡字符。
        extra = acb_mat([[frames[1][i,0],frames[2][i,0]] for i in range(19)])
        assert all(extra[i,j].is_zero() for i in untouched for j in (0,1))
        image = pivot*trivial+extra*acb_mat([[a,b],[0,d2.sqrt()]])
        for i in untouched:
            for j in (0,1):
                image[i,j]=weighted[i,j]
        minimum, margin = arb(1), arb(1)
        for row, w, q in zip(rows,weights,centers,strict=True):
            changed = row*image
            high = min((floor+sum(abs(z)**2 for z in r)).lower() for r in changed.tolist()[1:])
            assert high**2>p0
            error = (large(high**2,p0,p1)-reference).abs_upper()*sum(abs(v).upper() for v in w[1:])
            faces = [reference+w[0]*(weak(((q+sign*width)/mz)**2)-reference)+arb(0,error) for sign in (-1,1)]
            assert faces[0]>0>faces[1],(untouched,high,faces)
            minimum = min(minimum,high)
            margin = min(margin,faces[0].lower(),-faces[1].upper())
        print("保持不动的站位、其余截止下界、七对弱差面的最小余量",
              untouched,minimum.str(12),margin.str(12),flush=True)


if __name__=="__main__":
    root=Path(__file__).resolve().parent
    source=run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*source["verify"](root,return_operator="spectral"))
