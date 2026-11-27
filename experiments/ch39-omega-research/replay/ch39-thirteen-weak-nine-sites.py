"""九站共同酉同时消去第七至十九壳的十三个弱差。"""
from pathlib import Path
from runpy import run_path
from flint import arb, arb_mat, acb, acb_mat, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def compatible_rows(rows, radii):
    chosen = (0,1,3,6,7,9,11)
    other = (2,4,5,8,10,12)
    inverse = acb_mat([rows[j] for j in chosen]).inv()
    coefficients = acb_mat([[int(j == k) for j in range(13)] for k in other])
    predict = acb_mat([rows[j] for j in other])*inverse
    for h in range(6):
        for k,j in enumerate(chosen):
            coefficients[h,j] -= predict[h,k]
    angles = list(map(arb, ("-.4874122588629936","1.4185409399069038","-5.991740843796111",
        ".45285192963500354","-.770532230715546","-2.5793801095184263","-4.818273272016056",
        "3.7074649055911086","2.028202541790294","2.870254942621847","3.8900882122524925",
        "3.3491751537442025",".32327146508227367")))
    active = [j for j in range(13) if j != 4]
    width = arb(".00012")
    def system(expand):
        phases = [a+arb(0,width.upper()) if expand and j in active else a for j,a in enumerate(angles)]
        values = acb_mat([[r*acb(a.cos(),a.sin())] for r,a in zip(radii,phases,strict=True)])
        F = coefficients*values
        jac = [[acb(0,1)*coefficients[h,j]*values[j,0] for j in active] for h in range(6)]
        return (arb_mat([[getattr(z,part)] for part in ("real","imag") for z in F.entries()]),
                arb_mat([[getattr(z,part) for z in row] for part in ("real","imag") for row in jac]), values)
    F,J,_ = system(False)
    _,wide,values = system(True)
    preconditioner = J.mid().inv()
    identity = arb_mat([[int(i == j) for j in active] for i in active])
    contraction = max(sum(abs(z).upper() for z in row) for row in (identity-preconditioner*wide).tolist())
    displacement = max(abs(z).upper() for z in (preconditioner*F).entries())/width
    print("十二相位收缩常数、自映比", contraction.str(12), (contraction+displacement).str(12),flush=True)
    assert contraction+displacement < 1
    return inverse, acb_mat([[values[j,0]/10**8] for j in chosen])


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    graph = lib("three-ray-chiral-boundary")
    native = lib("native-rh-ray-pvm")
    frames,rotation = lib("cyclic-two-characters-nine")["characters"](graph["geo"]["sites"])
    lib("eleven-weak-eight-sites")["check_cyclic_geometry"](graph,rotation)
    frame = frames[0]
    gram = native["gram_cholesky"](root)[3]
    shell_rows, weights = [], []
    for shell in range(7,20):
        _,basis = lib("adjacent-pvm-transport")["spectrum"](graph,6,2_306_888_100_085+shell)
        shell_rows.append(adj(basis))
        weights.append([z.real for z in native["readout"](basis,range(19),[14,8,7,6],gram,3)])
    ctx.prec = 256
    rows = [(mat*frame).tolist()[0] for mat in shell_rows]
    assert all(sum(abs(z)**2 for z in row)>arb(".9") for row in rows)
    centers = list(map(arb,(*lib("six-cyclic-weak-horizon")["Q"],"506.579672393","504.414926074",
                       "510.694806713","520.456826617","525.342423153","535.119289427","535.106803223")))
    width = arb(".00005")
    floor, heavy, top = levels
    mp,mz = arb("1.22089e19"),arb("91.1876")
    radii = [10**8*((q+arb(0,width.upper()))/mp-floor).sqrt() for q in centers]
    inverse, selected = compatible_rows(rows,radii)
    means = (inverse*selected)*acb_mat([[1,0]])
    source = fifth*columns*acb_mat([[(heavy-floor).sqrt(),0],[0,(top-floor).sqrt()]])
    sites = (3,4,1,5,8,9,10,11,14)
    local = acb_mat([source.tolist()[j] for j in sites])
    average = adj(acb_mat([frame.tolist()[j] for j in sites]))
    dual = adj(average)*(average*adj(average)).inv()
    lift = dual*(means-adj(frame)*source+average*local)
    remaining = adj(local)*local - adj(lift)*lift
    a = remaining[0,0].real.sqrt()
    b = remaining[0,1]/a
    c2 = remaining[1,1].real-abs(b)**2
    assert remaining[0,0].real > arb(".01") and c2 > arb(".05")
    print("九站局部 Gram 的两个平方枢轴",remaining[0,0].real.str(12),c2.str(12),flush=True)
    extra = acb_mat([[(int(k==0)-int(k==1))/arb(2).sqrt(),(int(k==2)-int(k==3))/arb(2).sqrt()]
                     for k in range(9)])
    assert all(z.contains(0) for z in (average*extra).entries())
    image = acb_mat(source)
    for i,row in zip(sites,(lift+extra*acb_mat([[a,b],[0,c2.sqrt()]])).tolist(),strict=True):
        image[i,0],image[i,1] = row
    ref = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    low = native["weak_series"](56)
    large = lib("full-neutral-free-cutoff")["weak_column"]
    p0,p1 = (mz/mp)**2,(arb("1.176397e14")/mp)**2
    margins = []
    for mat,w,q in zip(shell_rows,weights,centers,strict=True):
        # 平凡字符的行直接从七条校准行恢复，使固定矩阵先相乘。
        trivial_rows = mat*frame
        projected = (trivial_rows*inverse)*selected
        lines = (mat*image).tolist()[1:]
        for k,row in enumerate(trivial_rows.tolist()[1:],start=1):
            if sum(abs(z)**2 for z in row)>arb(".5"):
                assert lines[k-1][0].overlaps(projected[k,0]) and lines[k-1][1].contains(0)
                lines[k-1] = [projected[k,0],acb(0)]
        assert all(z.is_finite() for row in lines for z in row)
        hi = min((floor+sum(max(arb(0),abs(z).lower())**2 for z in row)).lower() for row in lines)
        assert hi**2 > p0
        err = (large(hi**2,p0,p1)-ref).abs_upper()*sum(abs(z).upper() for z in w[1:])
        left = ref+w[0]*(low(((q-width)/mz)**2)-ref)
        right = ref+w[0]*(low(((q+width)/mz)**2)-ref)
        margins.extend(((left-err).lower(),(-right-err).lower()))
    assert all(m.is_finite() and m>arb("1e-8") for m in margins)
    print("十三对弱差面的最小余量",min(margins).str(12),flush=True)


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    original = run_path(str(root/"ch39-fifth-full-positive-joint.py"))["verify"]
    experiment(root,*original(root,return_operator="spectral"))
