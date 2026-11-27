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


def solve_matrix(root, data, centres, log_widths, *, return_transport=False, return_chart=False, tight_chart=False, component_box=False, name=NAME, upstream=None):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    sectors,incoming,readout,(low,h,k),(moments,heat,field),weights = data
    ctx.prec = 384
    identity = lib("rank-two-generator-up-quark")["eye"]
    eye = identity(19)
    coupling = lib("next-common-ckm")
    local = coupling["local"]
    seed = lib("one-generator-gravity-up-strange-bottom")["cayley"](root,identity,name+"-cayley.tsv")
    dirs = [(int(i),int(j),axis) for i,j,axis in (line.split() for line in
            (root/(name+"-directions.tsv")).read_text().splitlines())]
    assert len(set(dirs))==len(dirs)==26
    ov0 = readout*seed**3*incoming
    phases = [acb(ov0[i,0].real.mid(),ov0[i,0].imag.mid()) for i in LOW_ROWS]
    phases = [p/abs(p) for p in phases]
    zeros = [arb(0)]*26
    targets = list(map(arb,centres))
    wide_targets = [v*arb(0,arb(s).upper()).exp() for v,s in zip(targets,log_widths,strict=True)]

    def system(point,cutoffs, *, sources=None, source_tangents=None):
        sector_inputs, high_source = sources if sources is not None else (sectors,incoming)
        z,dz = local(seed,point,dirs,eye,source_tangents is None)
        z2=z*z
        d2=[e*z+z*e for e in dz]
        values,jac=[],[]
        for iq,((src,dst,w,cols,goal),op,ds) in enumerate(zip(sector_inputs,(adj(z2),z),([adj(d) for d in d2],dz),strict=True)):
            ov=dst*op*src
            output=arb_mat([[sum(v*abs(ov[i,j])**2 for j,v in enumerate(w))] for i in range(19)])
            values.extend((cols*output-goal).entries())
            derivatives=[]
            varied = ([dst*e*src for e in ds] if source_tangents is None else
                      [dst*op*de for de in source_tangents[iq]])
            for do in varied:
                dw=arb_mat([[2*sum(v*(ov[i,j].conjugate()*do[i,j]).real for j,v in enumerate(w))]
                            for i in range(19)])
                derivatives.append((cols*dw).entries())
            jac.extend(arb_mat(derivatives).transpose().tolist())

        ov=readout*z2*z*high_source
        d3=([readout*(e*z+z2*d)*high_source for e,d in zip(d2,dz,strict=True)] if source_tangents is None else
            [readout*z2*z*de for de in source_tangents[2]])
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
    projected=pre*f
    if upstream is not None:
        point_sources,tangents,widths=upstream
        f0,_,_=system(center,wide_targets,sources=point_sources)
        _,cross,_=system(center,wide_targets,source_tangents=tangents)
        # 先预条件化上游导数，再包围同一个上游 Cayley 参数在所有方程中的作用。
        projected=pre*f0+(pre*cross)*arb_mat([[arb(0,v.upper())] for v in widths])
    offset=max(abs(v).upper() for v in projected.entries())
    radius=2*offset.upper()
    radii=[radius]*26
    ident=arb_mat([[int(i==j) for j in range(26)] for i in range(26)])
    if component_box:
        residual=[abs(v).upper() for v in projected.entries()]
        radii=[2*max(v,radius/200) for v in residual]
        box=[v+arb(0,r.upper()) for v,r in zip(center,radii,strict=True)]
        _,wide,amps=system(box,wide_targets)
        defect=(ident-pre*wide).tolist()
        movement=[sum(abs(a).upper()*r for a,r in zip(row,radii,strict=True)) for row in defect]
        ratio=max((err+mov)/r for err,mov,r in zip(residual,movement,radii,strict=True))
        assert ratio<1
        print("逐坐标解盒的最大半宽与自映比上界",max(radii).str(13),ratio.upper().str(13),flush=True)
    else:
        box=[v+arb(0,radius.upper()) for v in center]
        _,wide,amps=system(box,wide_targets)
        lip=coupling["outer_sum"](ident-pre*wide)
        print("共同生成元中心偏移、参数盒半径、位移占比与导数误差行和",
              *[v.str(13) for v in (max(abs(a).upper() for a in center),radius,offset/radius,lip)],flush=True)
        assert lip+offset/radius<1
    sites = [a*a for a in amps]
    if return_transport or return_chart:
        z,_ = local(seed,center,dirs,eye,False)
        widths = [abs(b-v).upper() for b,v in zip(box,center,strict=True)]
        if return_chart:
            # 任一真根的各坐标偏移都受预条件中心误差和均值导数约束。
            # 迭代非负线性界仅收紧供下一壳使用的包围。
            if tight_chart:
                defect = ident-pre*wide
                magnitude = arb_mat([[abs(v).upper() for v in row] for row in defect.tolist()])
                residual = arb_mat([[abs(v).upper()] for v in projected.entries()])
                for _ in range(5):
                    bound = residual + magnitude*arb_mat([[v] for v in widths])
                    widths = [min(old, new.upper()) for old, new in zip(widths, bound.entries(), strict=True)]
                box = [v+arb(0,r.upper()) for v,r in zip(center,widths,strict=True)]
            wide,derivatives=local(seed,box,dirs,eye)
            return sites,(z,wide,derivatives,widths,sectors,incoming,(low,h,k))
        dz = 2*sum(a*a*(1 if i==j else 2) for a,(i,j,_) in zip(widths,dirs,strict=True)).sqrt()
        return sites,(z,dz,sectors,incoming,(low,h,k))
    return sites


def verify(root,*,return_transport=False,return_chart=False):
    adjust = run_path(str(root/"ch39-fourth-third-pivot-three-adjusters.py"))
    box = adjust["verify"](root)
    assert tuple(box["layers"])==LOW_ROWS
    result = solve_matrix(root,inputs(root),box["centers"],box["radii"],return_transport=return_transport,return_chart=return_chart)
    sites = result[0] if return_transport or return_chart else result
    high = [x for i,x in enumerate(sites) if i not in LOW_ROWS]
    shifts = [abs((x/c).log()).upper() for x,c in zip(high,box["high_centers"],strict=True)]
    assert max(shifts)<box["high_log_radius"]
    print("十六个高截止在联合盒内相对参照值的最大对数位移",max(shifts).str(14),flush=True)
    return result[1] if return_transport or return_chart else sites


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
