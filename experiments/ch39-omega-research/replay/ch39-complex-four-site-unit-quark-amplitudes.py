"""复四站共同味空间上，以两荷实际十九层严格介于零和一的振幅匹配同表中心。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


def adj(a):
    return a.conjugate().transpose()


def diag(values):
    return acb_mat([[x if i == j else 0 for j in range(len(values))] for i,x in enumerate(values)])


def unit(v):
    return v * (1/(adj(v)*v)[0,0].real.sqrt())


def four_basis():
    # 平凡块的复酉转角由两列高斯有理数的 Gram–Schmidt 唯一指定。
    first = unit(acb_mat([[acb(".98450852",".01345586")],[acb(".14041329",".10414428")]]))
    second = acb_mat([[acb(".16842528","-.04684922")],[acb("-.92382237","-.34057350")]])
    second = unit(second-first*(adj(first)*second)[0,0])
    rotate = acb_mat([[first[i,0],second[i,0]] for i in range(2)])
    h = 1/arb(2).sqrt()
    rows = [[rotate[i,0],rotate[i,1],0,0] for i in range(2)]
    return acb_mat(rows+[[0,0,h,h],[0,0,acb(0,h),acb(0,-h)]])


def coordinate_t(m):
    return [m[0,0].real,m[1,1].real,m[0,1].real,m[0,1].imag]


def reconstruct(rows, target, rotation, chosen):
    # 把实际四站谱投影写入所选平凡双重态／两个 C3 特征方向的基。
    native = [adj(rotation)*v*rotation for v in rows]
    assigned = [None]*19
    trivial = []
    classes = [[],[],[]]
    for i,p in enumerate(native):
        flags = [p[0,0].real+p[1,1].real>0, p[2,2].real>0, p[3,3].real>0]
        assert sum(flags)==1
        classes[flags.index(True)].append(i)
        if flags[0]:
            trivial.append(i)
            assigned[i]=arb("1e-9")
        else:
            assigned[i]=target[2 if flags[1] else 3,2 if flags[1] else 3].real
    assert [len(group) for group in classes]==[7,6,6] and all(i in trivial for i in chosen)
    coordinates=arb_mat([coordinate_t(native[i]) for i in chosen]).transpose()
    assert not coordinates.det().contains(0)
    current=sum((p*w for p,w in zip(native,assigned,strict=True)),acb_mat(4,4))
    rhs=arb_mat([[x] for x in coordinate_t(target-current)])
    increment=coordinates.solve(rhs,algorithm="precond")
    for i,v in zip(chosen,increment.entries(),strict=True):
        assigned[i]+=v
    assert all(0<w<1 for w in assigned)
    full=sum((p*w for p,w in zip(native,assigned,strict=True)),acb_mat(4,4))
    assert (full-target).contains(acb_mat(4,4))
    return full,assigned


def verify(root):
    fixed=run_path(str(root/"ch39-free-four-site-native-positive-spectrum.py"))
    ckm=run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    flavor=run_path(str(root/"ch39-three-gram-source-flavor.py"))
    sectors=fixed["native"](root)
    ctx.prec=320
    v=acb_mat(ckm["standard_matrix"](*[arb(i).sin() for i in (".22704",".04275",".00377")],arb("1.139")))
    up,down=[list(map(arb,s)) for s in ((".00000615",".00311",".8616"),(".0000135",".000268",".01401"))]
    a0=adj(v)*diag([v0/up[-1] for v0 in up])*v
    b0=diag([v0/down[-1] for v0 in down])
    # 原构造器给 z=i*残差。乘回 -i，使两部门的 yz 元同相且为实数。
    basis=fixed["target_basis"](a0,b0*arb(".1"))*diag([acb(0,1),1,acb(0,-1)])
    assert (adj(basis)*basis).contains(diag([1,1,1]))
    u,vv=arb(".00384"),arb(".0384")
    r=(vv/u).sqrt()
    a=(u*vv).sqrt()
    b=(1-a*a).sqrt()
    p=vv/(u+vv)
    f,g=p.sqrt(),(1-p).sqrt()
    embedding=acb_mat([[1,0,0],[0,0,b],[0,g,-a*f],[0,-f,-a*g]])
    rotation=four_basis()
    assert (adj(rotation)*rotation).contains(diag([1]*4))
    assert (adj(embedding)*embedding).contains(diag([1]*3))
    jproj=[]
    responses=[]
    for q,unscaled,masses,picks,rows in zip((-4,2),(a0,b0),(up,down),
                                             ((6,7,10,13),(6,7,10,12)),sectors,strict=True):
        m=adj(basis)*unscaled*basis*masses[-1]
        assert m[0,1].contains(0) and m[1,2].imag.contains(0)
        top,s,t = [m[j,j].real for j in range(3)]
        k,w=m[0,2],m[1,2].real
        char_plus,char_minus=s-w/u,s+w/vv
        block=acb_mat([[top,k/b],[k.conjugate()/b,(t-a*a*s-(u-vv)*w)/(b*b)]])
        eigen=sorted((v0.real for v0 in block.eig(algorithm="rump")),key=lambda v0:v0.mid())
        assert 0<eigen[0]<eigen[1]<1 and 0<char_plus<1 and 0<char_minus<1
        target=acb_mat([[block[0,0],block[0,1],0,0],[block[1,0],block[1,1],0,0],
                        [0,0,char_plus,0],[0,0,0,char_minus]])
        assert (adj(embedding)*target*embedding-m).contains(acb_mat(3,3))
        full,weights=reconstruct(rows,target,rotation,[i-1 for i in picks])
        lo,hi = (".9475540671066",".9475540671069") if q==-4 else (".3757768769212",".3757768769214")
        assert arb(lo)<max(w.lower() for w in weights) and max(w.upper() for w in weights)<arb(hi)
        responses.append(weights)
        compressed=adj(embedding)*full*embedding
        levels,proj=flavor["projections"](compressed*compressed)
        assert all(value.overlaps(mu**2) for value,mu in zip(levels,masses,strict=True))
        jproj.append(proj)
        print("荷、平凡块本征值、两个特征值",q,*[x.str(13) for x in eigen+[char_plus,char_minus]])
        print("19 个实际原壳物理振幅",*[x.str(13) for x in weights])
    pu,pd=jproj
    assert all((pu[i]*pd[j]).trace().real.overlaps(abs(v[i,j])**2) for i in range(3) for j in range(3))
    j=(pu[0]*pd[0]*pu[1]*pd[1]).trace().imag
    assert j>0 and j.overlaps((v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag)
    print("a、r、p、有符号 J",*[x.str(14) for x in (a,r,p,j)])
    return rotation*embedding, responses


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
