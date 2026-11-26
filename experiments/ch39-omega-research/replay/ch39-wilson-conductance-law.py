"""在原生联络各三角面上使用跨荷共用的非负 Wilson 权，检验保轻重塔的必要条件。"""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, acb, arb_mat, ctx

NUMBER = 2_306_888_100_085


def data(root):
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    recipe = run_path(str(root/"ch39-three-even-positive-towers.py"))
    ctx.prec=224
    height=acb.zeta_zero(NUMBER).imag
    n0=shell["laplacian"](0,height)
    native=shell["laplacian"](1,height)
    edges=[(i,j) for i,j,k,z in shell["geo"]["edges"]]
    base={e:abs(native[e]) for e in edges}
    phases={e:-native[e]/base[e] for e in edges}
    faces=[]
    for i,j,k in combinations(range(19),3):
        side=((i,j),(j,k),(i,k))
        if all(e in phases for e in side):
            loop=phases[i,j]*phases[j,k]*phases[i,k].conjugate()
            faces.append((side,loop))
    c=sum(n0[j,j].real for j in range(19))/18
    models={}
    for q in (2,4,6):
        L=shell["laplacian"](q,height)
        eigenvalues=L.eig(multiple=True)
        assert all(z.imag.contains(0) and z.real>0 for z in eigenvalues)
        eig=sorted([z.real for z in eigenvalues],key=float)
        assert all(a<b for a,b in zip(eig,eig[1:]))
        mass=sorted([v*abs(v-c) for v in eig],key=float)
        assert all(a<b for a,b in zip(mass,mass[1:])) and 4*mass[2]<c*c<128*mass[3]
        heavy=[(c+(c*c+128*v).sqrt())/2 for v in mass[3:]]
        lights=[[ (c-(c*c-4*v).sqrt())/2,(c+(c*c-4*v).sqrt())/2,(c+(c*c+4*v).sqrt())/2] for v in mass[:3]]
        opts=[(*sel,*heavy) for sel in product(*lights)]
        i,j,k=recipe["DATA"][q][2]
        target=[lights[i][0],lights[j][1],lights[k][2],*heavy]
        models[q]=(eig,mass,opts,target)
    return height, shell, c, edges, base, faces, models


def means(values):
    return sum(values),sum(v*v for v in values)


def verify(root):
    height,shell,c,edges,base,faces,models = data(root)
    assert len(edges)==135 and len(faces)==529
    # 整数面相位标签按实部辨识相同的三角面，包括反向的共轭。
    labels=run_path(str(root/"ch39-mincut-real-triangle-tomography.py"))["enumerate_faces"]()
    label=lambda side: labels[side[0][0],side[0][1],side[1][1]]
    index=69
    triple,loop=faces[index]
    assert triple==((1,8),(8,9),(1,9))
    smallest=4*loop.real**2
    same=[i for i,(side,phase) in enumerate(faces) if label(side)==label(triple)]
    assert same==[69,104,211]
    assert all(4*phase.real**2>smallest for i,(side,phase) in enumerate(faces) if i not in same)
    assert 1-(loop**2).real>0
    ratio=(smallest-1)**2
    assert smallest>arb("3.69031829") and ratio>arb("7.23781253")
    print("三角数、最小四荷/二荷增量比、最小六荷/二荷增量比、取等轨道",
          len(faces),smallest.str(17),ratio.str(17),same)
    baseline=2*sum(base.values())
    degrees=[sum(v for e,v in base.items() if i in e) for i in range(19)]
    second=sum(z*z for z in degrees)+2*sum(v*v for v in base.values())
    growth={}
    for q,(_,_,choices,selected) in models.items():
        first,Q=means(selected)
        choices=[means(a) for a in choices]
        # 根选择逐项独立，所有低支最小，所有高于中心的正支最大。
        lower=choices[0]
        upper=choices[-1]
        assert all(lower[0]<row[0]<upper[0] and lower[1]<row[1]<upper[1] for row in choices[1:-1])
        growth[q]=(first-baseline,lower[0]-baseline,upper[0]-baseline,
                   Q-second,lower[1]-second,upper[1]-second)
        print("荷、已选增迹、全轻根增迹范围、已选二阶增量、全轻根二阶增量范围",
              q,*[z.str(13) for z in growth[q]])
    for kind,two,six in (("指定原始根",growth[2][0],growth[6][0]),
                         ("允许全部轻根",growth[2][1],growth[6][2])):
        # 最优解只用一个极端面，两、六荷反向等误差。任一允许四荷目标均落在此误差内。
        error=(ratio*two-six)/(1+ratio)
        actual=(two+six)/(1+ratio)
        assert actual>0
        expected = "393.694846638581" if kind=="指定原始根" else "366.096125234543"
        assert abs(error-arb(expected))<arb("1e-11")
        assert all(abs(smallest*actual-w)<error for w in growth[4][1:3])
        coefficient=actual/(2*sum(base[e] for e in triple)*(1-(loop**2).real))
        print(kind,"三迹最大绝对误差的锐最小值、单面实现系数、逐特征值必要误差",
              error.str(17),coefficient.str(15),(error/19).str(15))
    # 二阶迹的增加是具有非负系数的线性项加二次项；顶点度及边权同时用柯西界。
    q2min=growth[2][4]
    q6max=growth[6][5]
    second_error=(ratio*q2min-q6max)/(1+ratio)
    corrected=ratio*q2min+ratio*(ratio-1)*(arb(1)/19+arb(1)/270)*growth[2][1]**2
    assert second_error>arb("24574.62") and corrected-q6max>arb("816665.53")
    print("仅用二阶迹的统一误差下界、若二荷两矩达标则六荷二阶迹至少超目标",
          second_error.str(17),(corrected-q6max).str(17))
    exponent=(growth[6][2]/growth[2][1]).log()/ratio.log()
    assert exponent<arb(".049")
    print("任意轻根时共同幂响应的指数必要上界",exponent.str(12))
    # 对已发表的荷六分支，目标最低根就是原最低根，且原最低特征向量在每个顶点非零。
    matrix=shell["laplacian"](6,height)
    values,vectors=matrix.eig(right=True,algorithm="rump")
    j=min(range(19),key=lambda i:float(values[i].real))
    assert all(values[j].real<z.real for i,z in enumerate(values) if i!=j)
    lam=values[j].real
    assert 0<lam<c/2
    folded=[(i,z.real*abs(z.real-c)) for i,z in enumerate(values)]
    ordered=sorted(folded,key=lambda v:float(v[1]))
    assert all(a[1]<b[1] for a,b in zip(ordered,ordered[1:])) and ordered[2][0]==j
    assert (models[6][3][0]-lam).contains(0)
    norm=sum(abs(vectors[i,j])**2 for i in range(19)).sqrt()
    ground_min=min((abs(vectors[i,j])/norm).lower() for i in range(19))
    assert ground_min>arb(".22024")
    print("原荷六最低根、归一基态最小顶点振幅、保荷二目标时荷一必要增迹",
          lam.str(17),ground_min.str(14),(growth[2][1]/4).str(14))


if __name__=="__main__":
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
