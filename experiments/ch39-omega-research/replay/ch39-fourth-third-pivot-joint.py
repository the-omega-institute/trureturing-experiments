"""以共同第三步与共同三代平面同时配平夸克读数及第三枢轴的荷六方程。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


NAME = "ch39-fourth-third-pivot-joint"
LOW_ROWS = (9, 0, 16)


def adj(matrix):
    return matrix.conjugate().transpose()


def inputs(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    sectors,incoming,readout,levels,thermal = lib("fourth-full-ckm-gravity")["verify"](root,return_inputs=True)
    _,_,pivots,_ = next(row for row in lib("native-rh-ray-pvm")["shell_weights"](root,2_306_888_100_089) if row[0]==6)
    return sectors,incoming,readout,levels,thermal,[q.real for q in pivots[2]]


def solve_matrix(root, data, centres, log_widths):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    sectors,incoming,readout,(low,h,k),(moments,heat,field),weights = data
    ctx.prec = 384
    identity = lib("rank-two-generator-up-quark")["eye"]
    eye = identity(19)
    coupling = lib("next-common-ckm")
    local = coupling["local"]
    seed = lib("one-generator-gravity-up-strange-bottom")["cayley"](root,identity,NAME+"-cayley.tsv")
    dirs = [(int(i),int(j),axis) for i,j,axis in (line.split() for line in
            (root/(NAME+"-directions.tsv")).read_text().splitlines())]
    assert len(set(dirs))==len(dirs)==26
    ov0 = readout*seed**3*incoming
    phases = [acb(ov0[i,0].real.mid(),ov0[i,0].imag.mid()) for i in LOW_ROWS]
    phases = [p/abs(p) for p in phases]
    zeros = [arb(0)]*26
    targets = list(map(arb,centres))
    wide_targets = [v*arb(0,arb(s).upper()).exp() for v,s in zip(targets,log_widths,strict=True)]

    def system(point,cutoffs):
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
        amps=[low+(h-low)*abs(ov[i,0])**2+(k-low)*abs(ov[i,1])**2 for i in range(19)]
        for i,x,phase in zip(LOW_ROWS,cutoffs,phases,strict=True):
            amps[i]=x.sqrt()
            eta=((amps[i]-low)/(h-low)).sqrt()
            res=ov[i,0]-eta*phase
            values.extend((res.real,res.imag,ov[i,1].real,ov[i,1].imag))
            jac.extend([[getattr(d[i,j],part) for d in d3] for j in (0,1) for part in ("real","imag")])

        thermal_values=[-arb(1),-arb("1e-122")]
        slopes=[[],[]]
        for i,a in enumerate(amps):
            for j,v in enumerate(moments(field,heat,a*a)):
                thermal_values[j]+=weights[i]*v
            if i in LOW_ROWS:
                for col in slopes:
                    col.append(arb(0))
                continue
            assert a>0
            ex=[(-mass/a**2).exp() for _,_,mass in field]
            slopes[0].append(weights[i]*a*sum(c*v for (_,c,_),v in zip(field,ex,strict=True))/arb.pi())
            slopes[1].append(-weights[i]*a**3*sum(g*v for (g,_,_),v in zip(field,ex,strict=True))/(2*arb.pi()))
        derivative=arb_mat([[2*((h-low)*(ov[i,0].conjugate()*d[i,0]).real+
                                (k-low)*(ov[i,1].conjugate()*d[i,1]).real) for d in d3] for i in range(19)])
        values.extend(thermal_values)
        jac.extend((arb_mat(slopes)*derivative).tolist())
        return arb_mat([[v] for v in values]),arb_mat(jac),amps

    center=zeros
    for _ in range(3):
        f,j,_=system(center,targets)
        delta=j.mid().solve(f.mid()).entries()
        center=[(v-d).mid() for v,d in zip(center,delta,strict=True)]
    f,j,_=system(center,wide_targets)
    pre=j.mid().inv()
    assert not pre.det().contains(0)
    offset=max(abs(v).upper() for v in (pre*f).entries())
    radius=2*offset.upper()
    box=[v+arb(0,radius.upper()) for v in center]
    _,wide,amps=system(box,wide_targets)
    ident=arb_mat([[int(i==j) for j in range(26)] for i in range(26)])
    lip=coupling["outer_sum"](ident-pre*wide)
    print("共同生成元中心偏移、参数盒半径、位移占比与导数误差行和",
          *[v.str(13) for v in (max(abs(a).upper() for a in center),radius,offset/radius,lip)],flush=True)
    assert lip+offset/radius<1
    return [a*a for a in amps]


def verify(root):
    adjust = run_path(str(root/"ch39-fourth-third-pivot-three-adjusters.py"))
    box = adjust["verify"](root)
    assert tuple(box["layers"])==LOW_ROWS
    sites = solve_matrix(root,inputs(root),box["centers"],box["radii"])
    high = [x for i,x in enumerate(sites) if i not in LOW_ROWS]
    shifts = [abs((x/c).log()).upper() for x,c in zip(high,box["high_centers"],strict=True)]
    assert max(shifts)<box["high_log_radius"]
    print("十六个高截止在联合盒内相对参照值的最大对数位移",max(shifts).str(14),flush=True)
    return sites


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
