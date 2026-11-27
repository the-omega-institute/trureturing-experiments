"""在前两步的整个联合盒上用一个第三步生成元与共同三平面实现 N+4 的全部目标。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


QUATERNION = (-12714083761,-70309030904,-460303187067,884881895070)


def adj(a):
    return a.conjugate().transpose()


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    x,y,incoming,readout,(low,h,k),weights,thermal = lib("fourth-shell-two-gravity-rules")["shell_data"](root)
    moments,heat,field = thermal
    original = lib("one-generator-gravity-full-ckm")
    emb = lib("common-embedding-central-right-ckm")
    independent = lib("independent-unitary-quark-right-holdout")
    rank = lib("rank-two-generator-up-quark")
    n = rank["N"]
    _,training,_ = emb["verify_nodes"](root,((n+1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    targets = original["targets"](root,independent,emb)
    make_block = lib("next-coherent-full-ckm")["block_from_coordinates"]
    coords = independent["coordinates"]
    rows = lib("free-four-site-native-positive-spectrum")["native"](root,n+4)
    ctx.prec = 384
    length = arb(sum(v*v for v in QUATERNION)).sqrt()
    qx,qy,qz,qw = [arb(v)/length for v in QUATERNION]
    rotation = acb_mat([[acb(qw,-qz),acb(-qy,-qx),0,0],[acb(qy,-qx),acb(qw,qz),0,0],
                        [0,0,1,0],[0,0,0,1]])
    _,_,f,r0 = emb["common_embedding"](*independent["GEOMETRY"])
    frame = rotation*r0*f
    assert (adj(frame)*frame).contains(rank["eye"](3))
    spec = lib("adjacent-pvm-transport")["spectrum"]
    shell = lib("three-ray-chiral-boundary")
    sectors = []
    for q,old,projectors,propagation in zip((-4,2),targets,rows,(adj(y*y)*adj(x*x),y*x),strict=True):
        block = make_block(old.entries())
        goal = rotation*block*adj(rotation)
        assert (adj(frame)*goal*frame-adj(r0*f)*block*r0*f).contains(acb_mat(3,3))
        src,dst = [spec(shell,q,j)[1] for j in (n+1,n+4)]
        sectors.append((propagation*src,adj(dst),training[n+1,q],arb_mat([coords(p) for p in projectors]).transpose(),
                        arb_mat([[v] for v in coords(goal)])))
    coupling = lib("next-common-ckm")
    local = coupling["local"]
    eye = rank["eye"](19)
    seed = lib("one-generator-gravity-up-strange-bottom")["cayley"](
        root,rank["eye"],"ch39-fourth-full-ckm-gravity-cayley.tsv")
    dirs = []
    for line in (root/"ch39-fourth-full-ckm-gravity-directions.tsv").read_text().splitlines():
        i,j,axis = line.split()
        dirs.append((int(i),int(j),axis))
    assert len(set(dirs))==len(dirs)==18
    link = (readout*seed**3*incoming)[15,0]
    # 沿中点的固定相位选取该复重叠的代表。
    link = acb(link.real.mid(),link.imag.mid())
    phase = link/abs(link)
    native = lib("native-rh-ray-pvm")
    weak = native["weak_series"](56)
    high_weak = lib("full-neutral-free-cutoff")["weak_column"]
    mz,mp = arb("91.1876"),arb("1.22089e19")
    p0,p1 = (mz/mp)**2,(arb("1.176397e14")/mp)**2
    left,right = arb("1.4"),arb("2.1")
    def rough(t):
        return weights[15]*weak(t)+(1-weights[15])*arb("-.074256577194663365545")
    assert rough(left)>0>rough(right)
    for _ in range(44):
        mid=(left+right)/2
        if rough(mid)>0:
            left=mid
        else:
            assert rough(mid)<0
            right=mid
    left-=arb("1e-11")
    right+=arb("1e-11")
    time=left.union(right)

    def system(point,t):
        z,dz = local(seed,point,dirs,eye)
        z2=z*z
        d2=[e*z+z*e for e in dz]
        values,jac=[],[]
        for (src,dst,w,cols,goal),op,ds in zip(sectors,(adj(z2),z),([adj(d) for d in d2],dz),strict=True):
            ov=dst*op*src
            output=arb_mat([[sum(v*abs(ov[i,j])**2 for j,v in enumerate(w))] for i in range(19)])
            values.extend((cols*output-goal).entries())
            derivatives=[]
            for e in ds:
                do=dst*e*src
                dw=arb_mat([[2*sum(v*(ov[i,j].conjugate()*do[i,j]).real for j,v in enumerate(w))]
                            for i in range(19)])
                derivatives.append((cols*dw).entries())
            jac.extend(arb_mat(derivatives).transpose().tolist())
        ov=readout*z2*z*incoming
        d3=[readout*(e*z+z2*d)*incoming for e,d in zip(d2,dz,strict=True)]
        cut=mz*t.sqrt()/mp
        eta=((cut-low)/(h-low)).sqrt()
        res=ov[15,0]-eta*phase
        values.extend((res.real,res.imag,ov[15,1].real,ov[15,1].imag))
        jac.extend([[getattr(d[15,j],part) for d in d3] for j in (0,1) for part in ("real","imag")])
        amps=[low+(h-low)*abs(ov[i,0])**2+(k-low)*abs(ov[i,1])**2 for i in range(19)]
        amps[15]=cut
        # 在满足四个重叠方程时，第16层与热核计算中的 cut 完全相同。
        terms=list(moments(field,heat,cut*cut))
        thermal_values=[weights[15]*a-target for a,target in zip(terms,(1,arb("1e-122")),strict=True)]
        slopes=[[],[]]
        for i,a in enumerate(amps):
            if i==15:
                for col in slopes:
                    col.append(arb(0))
                continue
            assert a>arb("1e-7")
            for j,v in enumerate(moments(field,heat,a*a)):
                thermal_values[j]+=weights[i]*v
            ex=[(-mass/a**2).exp() for _,_,mass in field]
            slopes[0].append(weights[i]*a*sum(c*v for (_,c,_),v in zip(field,ex,strict=True))/arb.pi())
            slopes[1].append(-weights[i]*a**3*sum(g*v for (g,_,_),v in zip(field,ex,strict=True))/(2*arb.pi()))
        derivative=arb_mat([[2*((h-low)*(ov[i,0].conjugate()*d[i,0]).real+
                                (k-low)*(ov[i,1].conjugate()*d[i,1]).real) for d in d3] for i in range(19)])
        values.extend(thermal_values)
        jac.extend((arb_mat(slopes)*derivative).tolist())
        return arb_mat([[v] for v in values]),arb_mat(jac),amps

    center=[arb(0)]*18
    for _ in range(3):
        f,j,_=system(center,time.mid())
        delta=j.mid().solve(f.mid())
        center=[(z-d).mid() for z,d in zip(center,delta.entries(),strict=True)]
    f,j,_=system(center,time)
    pre=j.mid().inv()
    assert not pre.det().contains(0)
    offset=max(abs(v).upper() for v in (pre*f).entries())
    radius=2*offset.upper()
    box=[v+arb(0,radius.upper()) for v in center]
    _,wide,amps=system(box,time)
    ident=arb_mat([[int(i==j) for j in range(18)] for i in range(18)])
    lip=coupling["outer_sum"](ident-pre*wide)
    print("第三步最大中心、半宽、中心位移占比、导数行和",*[v.str(12) for v in
          (max(abs(a).upper() for a in center),radius,offset/radius,lip)],flush=True)
    assert lip+offset/radius<1
    other=sum(c*high_weak(a*a,p0,p1) for i,(c,a) in enumerate(zip(weights,amps,strict=True)) if i!=15)
    ends=[weights[15]*weak(t)+other for t in (left,right)]
    assert ends[0]>0>ends[1]
    print("第16层截止 GeV、弱差两端、另外十八层最小振幅",
          *[v.str(13) for v in (mz*time.sqrt(),*ends,min(a.lower() for i,a in enumerate(amps) if i!=15))],
          flush=True)
    return center,radius,time,frame


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
