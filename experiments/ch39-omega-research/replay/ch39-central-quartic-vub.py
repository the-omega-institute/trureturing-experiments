"""中央实际壳以一条有正底的四次响应在同一空间实现六中心质量与 |Vub|。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, arb_series, ctx


def rooted_shape(model,width,guesses,nodes,rows,masses):
    coeff = [arb(x) for x in (".9979171475059213",".002839084023218489",
                              "19.13848104749016","8.056457264569582")]+[width]
    center = list(guesses)
    for _ in range(5):
        f,j = model["constraints"](coeff,center,nodes,rows,masses)
        step = j.solve(f)
        center = [(x-y).mid() for x,y in zip(center,step.entries(),strict=True)]
    radius = arb("1e-35")
    box = [x+arb(0,radius.upper()) for x in center]
    f,j = model["constraints"](coeff,center,nodes,rows,masses)
    inv = j.inv()
    _,jac = model["constraints"](coeff,box,nodes,rows,masses)
    err = arb_mat([[1,0],[0,1]])-inv*jac
    gain = max(sum(abs(x).upper() for x in row) for row in err.tolist())
    assert max(abs(x).upper() for x in (inv*f).entries())+radius*gain < radius
    residue,c = model["read"](coeff,box,nodes,rows,masses)
    assert all(0<x<1 for row in residue for x in row)
    sx,sy = [(rr[1]-rr[0])/(rr[1]+rr[0]) for rr in residue]
    assert (1-c*c)*(1-sx*sx)-(sy-c*sx)**2 > 0
    return box,model["ckm_branches"](coeff,box,nodes,rows,masses,verbose=False)


def verify(root):
    native = run_path(str(root/"ch39-free-four-site-native-positive-spectrum.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ordered = run_path(str(root/"ch39-bounded-common-rational-amplitude.py"))["ordered"]
    model = run_path(str(root/"ch39-common-quartic-six-masses.py"))
    n = 2_306_888_100_085
    ps = native["native"](root,n)
    ctx.prec,ctx.cap = 320,2
    masses = [list(map(arb,k)) for k in ((".00000615",".00311",".8616"),(".0000135",".000268",".01401"))]
    height = acb.zeta_zero(n).imag
    nodes = [ordered(shell,q,height) for q in (-4,2)]
    rows = []
    for sector in ps:
        coord = [(p[0,0].real,p[1,1].real,p[0,1].real,p[0,1].imag,
                  p[2,2].real+p[2,3].imag,p[2,2].real-p[2,3].imag) for p in sector]
        rows.append(list(zip(*coord,strict=True)))
    base = arb(".19818774905676692")
    guesses = [arb("19.088219245926343"),arb(".24650474619647347")]
    wanted = arb(".00377").sin()
    # 割线迭代仅提供候选中心。存在性由最终两个严格符号和统一压缩给出。
    widths = [base,base+arb("1e-9")]
    errors = []
    for _ in range(5):
        widths = widths[-2:]
        errors = []
        for w in widths:
            _,response = rooted_shape(model,w,guesses,nodes,rows,masses)
            errors.append((response[0][0]-wanted).mid())
        if abs(errors[1])<arb("1e-25"):
            break
        widths.append((widths[1]-errors[1]*(widths[1]-widths[0])/(errors[1]-errors[0])).mid())
    wc=widths[-1]
    interval = wc+arb(0,arb("1e-17").upper())
    exact_mid,_ = rooted_shape(model,wc,guesses,nodes,rows,masses)
    center = [x.mid() for x in exact_mid]
    radii = [arb("1e-17"),arb("1e-13")]
    box = [x+arb(0,r.upper()) for x,r in zip(center,radii,strict=True)]
    constants=[arb(x) for x in (".9979171475059213",".002839084023218489",
                                "19.13848104749016","8.056457264569582")]
    assert interval>0 and all(not (h-constants[3]).contains(0) for h in box)
    f,j = model["constraints"](constants+[wc],center,nodes,rows,masses)
    inv=j.inv()
    rr,_=model["read"](constants+[arb_series([interval,1])],center,nodes,rows,masses)
    variation=inv*arb_mat([[(a-b)[1]] for a,b in zip(rr[0][2:],rr[1][2:],strict=True)])
    _,jj=model["constraints"](constants+[interval],box,nodes,rows,masses)
    defect=arb_mat([[1,0],[0,1]])-inv*jj
    ratios=[]
    for i,row in enumerate(defect.tolist()):
        shift=abs((inv*f)[i,0])+abs(variation[i,0])*(interval-wc).abs_upper()
        ratios.append((shift+sum(abs(a)*q for a,q in zip(row,radii,strict=True)))/radii[i])
    assert max(z.upper() for z in ratios)<1
    # 参数化的不动点在整个区间连续。两端所用单点根必须属于同一个盒。
    ends=[]
    for w in (interval.lower(),interval.upper()):
        point,response=rooted_shape(model,w,center,nodes,rows,masses)
        assert all(B.contains(v) for B,v in zip(box,point,strict=True))
        ends.append(response[0][0]-wanted)
    assert ends[0]>0>ends[1]
    # 上述区间内两条法向分支均存在且无相切。
    rr,c=model["read"](constants+[interval],box,nodes,rows,masses)
    assert all(0<x<1 for row in rr for x in row)
    su,sd=[(q[1]-q[0])/(q[1]+q[0]) for q in rr]
    assert (1-c*c)*(1-su*su)-(sd-c*su)**2>0
    branches=model["ckm_branches"](constants+[interval],box,nodes,rows,masses,verbose=False)
    for vub,j,p in branches:
        vcb=p[1][2].sqrt()
        assert vcb>arb(".998") and abs(j)<arb("1.5e-8")
        print("参数全盒的 |Vcb| 和有符号 J",vcb.str(11),j.str(10),flush=True)
    print("形状参数、峰位、统一压缩界、两端 |Vub| 残差",interval.str(17),
          *[x.str(15) for x in box],max(v.upper() for v in ratios).str(10),
          *[v.str(9) for v in ends],flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
