"""固定同一个八站酉同时消去第七至十七壳弱差。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def check_cyclic_geometry(graph, rotation):
    """以整数点积和三重积核对全部有向边在三循环下的几何数据。"""
    geo = graph["geo"]
    vertices = geo["sites"]
    pole_index = graph["soft"]["pole"]
    pole = vertices[pole_index]
    assert pole_index == rotation[pole_index] == 14
    assert pole == (17, 17, 17)
    tuned = set(graph["soft"]["star"]["orbit"])
    assert {rotation[i] for i in tuned} == tuned
    lengths = {(i, j): dist for i, j, dist, _ in geo["edges"]}
    assert len(lengths) == 135
    for (i, j), dist in lengths.items():
        i1, j1 = rotation[i], rotation[j]
        assert lengths[tuple(sorted((i1, j1)))] == dist
        # 原输入按从小到大编号附加时钟，仅两个带钟距离需要检验朝向。
        assert dist not in (724, 804) or i1 < j1
        x, y, rx, ry = (vertices[k] for k in (i, j, i1, j1))
        dot, det = geo["dot"], geo["det"]
        assert dot(pole, x) == dot(pole, rx) and dot(y, pole) == dot(ry, pole)
        assert dot(x, y) == dot(rx, ry) and det(pole, x, y) == det(pole, rx, ry)


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    ctx.prec=256
    graph=lib("three-ray-chiral-boundary")
    native=lib("native-rh-ray-pvm")
    frames, rotation = lib("cyclic-two-characters-nine")["characters"](graph["geo"]["sites"])
    frame = frames[0]
    check_cyclic_geometry(graph, rotation)
    weights_gram=native["gram_cholesky"](root)[3]
    floor,middle,top=levels
    source=fifth*columns*acb_mat([[(middle-floor).sqrt(),0],[0,(top-floor).sqrt()]])
    projected=adj(frame)*source
    pole=acb_mat([projected.tolist()[6]])
    pole_norm=sum(abs(z)**2 for z in pole.entries()).sqrt()
    assert pole_norm > arb(".02")
    direction=pole/pole_norm
    rows=[]
    weights=[]
    for j in range(7,18):
        _,basis=lib("adjacent-pvm-transport")["spectrum"](graph,6,2_306_888_100_085+j)
        rows.append(adj(basis))
        weights.append([z.real for z in native["readout"](basis,range(19),[14,8,7,6],weights_gram,3)])
    b=[(row*frame).tolist()[0] for row in rows]
    # 简单谱与三循环对易，每条本征线只有一个字符。严格正的平凡分量确定其字符。
    assert all(sum(max(arb(0),abs(z).lower())**2 for z in row) > arb(".9") for row in b)
    assert all(w[0] > 0 for w in weights)
    chosen=[1,3,4,6,8,10]
    extra=[0,2,5,7,9]
    inv=acb_mat([b[j][:6] for j in chosen]).inv()
    synth=acb_mat([b[j][:6] for j in extra])*inv
    C=acb_mat([[int(j==ex) for j in range(11)] for ex in extra])
    for i in range(len(extra)):
        for k,j in enumerate(chosen):
            C[i,j]-=synth[i,k]
    scale=arb(10)**8
    d=(acb_mat([[b[j][6]] for j in extra])-synth*acb_mat([[b[j][6]] for j in chosen]))*(pole_norm*scale)
    centers=list(map(arb,(*lib("six-cyclic-weak-horizon")["Q"],"506.579672393","504.414926074","510.694806713","520.456826617","525.342423153")))
    width=arb(".00005")
    planck,mz=arb("1.22089e19"),arb("91.1876")
    r=[((q+arb(0,width.upper()))/planck-floor).sqrt()*scale for q in centers]
    angles=list(map(arb,('-2.8334627213829195', '-4.522787132657973', '3.5281940505351215', '2.849927995182405', '1.8082708731402373', '-2.0479103405920087', '-3.110986559454634', '3.0041358500245106', '-1.9340893977735696', '-2.246885664597755', '-1.0355624123459037')))
    active=[0,1,2,3,4,5,7,8,9,10]
    phase_radius=arb(".0001")
    def system(wide):
        phase=[t+arb(0,phase_radius.upper()) if wide and k in active else t for k,t in enumerate(angles)]
        a=acb_mat([[rk*acb(t.cos(),t.sin())] for rk,t in zip(r,phase,strict=True)])
        defect=C*a-d
        df=[[C[k,j]*acb(0,1)*a[j,0] for j in active] for k in range(len(extra))]
        J=arb_mat([[getattr(z,part) for z in row] for part in ('real','imag') for row in df])
        F=arb_mat([[getattr(z,part)] for part in ('real','imag') for z in defect.entries()])
        return F,J,a
    value,jmid,_=system(False)
    _,jac,targets=system(True)
    pre=arb_mat([[z.mid() for z in row] for row in jmid.tolist()]).inv()
    I=arb_mat([[int(i==j) for j in range(len(active))] for i in range(len(active))])
    lip=max(sum(abs(z).upper() for z in row) for row in (I-pre*jac).tolist())
    move=max(abs(z).upper() for z in (pre*value).entries())
    ratio=move/phase_radius+lip
    print("十相位压缩常数、自映比",lip.str(10),ratio.str(10),flush=True)
    assert lip.is_finite() and ratio.is_finite() and lip < arb(".0022") and ratio < arb(".07")
    # 对每个截止盒中的实向量，压缩映射给出唯一且连续的十个相位。
    # 在解上，选定六条谱行决定六个自由均值；余下五行满足复相容方程。
    first6=(inv*acb_mat([[targets[j,0]/scale] for j in chosen])
            -(inv*acb_mat([[b[j][6]] for j in chosen]))*pole_norm)
    avg_trivial=first6*direction
    sites=(1,3,4,5,8,9,10,11)
    groups=[[i for i in sites if abs(frame[i,k])>arb(".5")] for k in range(6)]
    assert tuple(map(len,groups))==(2,2,1,1,1,1)
    assert set(sites)==set().union(*groups) and 14 not in sites
    local=acb_mat([source.tolist()[i] for i in sites])
    avg=acb_mat([[sum(source[i,h] for i in g)/arb(len(g)).sqrt()
                  +(arb(3)/len(g)).sqrt()*(avg_trivial[k,h]-projected[k,h])
                  for h in (0,1)] for k,g in enumerate(groups)])
    rem=adj(local)*local-adj(avg)*avg
    a=rem[0,0].real.sqrt()
    b12=rem[0,1]/a
    d2=rem[1,1].real-abs(b12)**2
    print("局部 Gram 两个 Cholesky 平方枢轴",rem[0,0].real.str(10),d2.str(10),flush=True)
    assert rem[0,0].real>arb(".0214") and d2>arb(".1318")
    comp=acb_mat([[a,b12],[0,d2.sqrt()]])
    image=acb_mat(source)
    kextra=0
    for k,g in enumerate(groups):
        for i in g:
            for h in (0,1):
                image[i,h]=avg[k,h]/arb(len(g)).sqrt()
        if len(g)==2:
            for h in (0,1):
                v=comp[kextra,h]/arb(2).sqrt()
                image[g[0],h]+=v
                image[g[1],h]-=v
            kextra+=1
    assert kextra==2
    weak=native["weak_series"](56)
    large=lib("full-neutral-free-cutoff")["weak_column"]
    ref=arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    p0,p1=(mz/planck)**2,(arb("1.176397e14")/planck)**2
    margins=[]
    for row,w,q in zip(rows,weights,centers,strict=True):
        high=min((floor+sum(max(arb(0),abs(z).lower())**2 for z in line)).lower()
                 for line in (row*image).tolist()[1:])
        assert high.is_finite() and high**2>p0
        err=(large(high**2,p0,p1)-ref).abs_upper()*sum(abs(z).upper() for z in w[1:])
        m1=(ref+w[0]*(weak(((q-width)/mz)**2)-ref)-err).lower()
        m2=(-ref-w[0]*(weak(((q+width)/mz)**2)-ref)-err).lower()
        assert m1.is_finite() and m2.is_finite() and min(m1,m2)>0,(len(margins)+7,high,err,m1,m2)
        margins.append(min(m1,m2))
    # 以上为十一维截止盒的相对两面的严格异号。Poincaré–Miranda 给出共同零点。
    assert min(margins)>arb("1.4e-8")
    print("第七至十七壳十一对弱差面最小余量",min(margins).str(12),flush=True)
    return ratio, min(margins)


if __name__=="__main__":
    root=Path(__file__).resolve().parent
    source=run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*source["verify"](root,return_operator="spectral"))
