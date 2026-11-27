"""检验指定三壳弱截止时，移动第五壳重列最少的唯一标架。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def adj(a):
    return a.conjugate().transpose()


def positive_root(a, identity):
    # 正定二阶矩阵满足 Cayley–Hamilton，故正平方根只需迹和行列式。
    d = a.det().real
    tr = (a[0,0]+a[1,1]).real
    assert a[0,0].real > 0 and d > 0
    return (a + identity*d.sqrt())*(1/(tr+2*d.sqrt()).sqrt())


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    source = fifth*columns
    shells = []
    for n in range(2_306_888_100_092,2_306_888_100_095):
        _,_,_,base,_,weights,heat,field = lib("fifth-seventh-shortest-site-gates")["geometry"](
            root,fifth,columns,levels,number=n)
        shells.append((adj(base),weights,heat,field))
    ctx.prec = 256
    rows = acb_mat([[left[15,i] for i in range(19)] for left,_,_,_ in shells])
    eye = acb_mat([[1,0],[0,1]])
    direct = rows*source
    gram = rows*adj(rows)
    assert not gram.det().contains(0)
    dual = adj(rows)*gram.inv()
    projector = acb_mat([[int(i==j) for j in range(19)] for i in range(19)])-dual*rows
    projected = projector*source
    low = levels[0]
    planck,mz = arb("1.22089e19"),arb("91.1876")
    q0 = list(map(arb,lib("common-three-shell-nine-endpoints")["Q_CENTERS"]))
    width = arb(".000001")
    q = [v+arb(0,width) for v in q0]
    scales = [((v/planck-low)/sum((levels[j+1]-low)*abs(direct[i,j])**2 for j in range(2))).sqrt()
              for i,v in enumerate(q)]
    assert all(0 < s < 1 for s in scales)
    prescribed = acb_mat([[scales[i]*direct[i,j] for j in range(2)] for i in range(3)])
    lift = dual*prescribed
    t = positive_root(eye-adj(lift)*lift,eye)
    z = projected*t
    cov = positive_root(adj(z)*z,eye)
    correction = t*cov.inv()*t
    weak = lib("native-rh-ray-pvm")["weak_series"](56)
    large = lib("full-neutral-free-cutoff")["weak_column"]
    ref = arb("-.074256577194663365545")
    p0,p1 = (mz/planck)**2,(arb("1.176397e14")/planck)**2
    targets = (arb(1),arb("1e-122"))
    moments = lib("signed-cutoff-gravity")["moments"]
    for i,(left,w,heat,field) in enumerate(shells):
        p = (left*dual)*prescribed+((left*projector)*source)*correction
        amplitudes = [low+sum((levels[j+1]-low)*abs(p[k,j])**2 for j in range(2)) for k in range(19)]
        remaining = [(w[k],a.lower()) for k,a in enumerate(amplitudes) if k!=15]
        minimum = min(a for _,a in remaining)
        assert minimum > 0 and minimum**2 > p0 and sum(w).contains(1)
        error = sum(weight*(large(a**2,p0,p1)-ref) for weight,a in remaining)
        faces = [ref+w[15]*(weak(((q0[i]+sgn*width)/mz)**2)-ref)+error for sgn in (-1,1)]
        print("壳、其余振幅下界、两面弱差",i+7,minimum.str(12),faces[0].str(12),faces[1].str(12),flush=True)
        assert faces[0] > 0 > faces[1]
        # 第 16 线的振幅按输入精确施加。
        amplitudes[15] = q[i]/planck
        _,column = lib("fifth-seventh-two-pair-thermal-rotation")["thermal_columns"](field,moments,heat)
        readings = [column(a) for a in amplitudes]
        values = [sum(weight*r[k] for weight,r in zip(w,readings,strict=True)) for k in range(2)]
        print("牛顿、真空读数端点",[(v.lower().str(14),v.upper().str(14)) for v in values],flush=True)
        assert all(not v.overlaps(t) for v,t in zip(values,targets,strict=True))


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    original = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    experiment(root,*original["verify"](root,return_operator="spectral"))
